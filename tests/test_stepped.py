import ctypes
import json
import math

import numpy as np
import pytest
from build123d import Cylinder, GeomType, Pos, import_step, import_stl
from OCP.BRepExtrema import BRepExtrema_DistShapeShape
from test_planetary import GAP_MAX, MARGIN, _faces

from reducer import BEARINGS, MIN, interference
from stepped import (
    Params,
    build,
    design_warnings,
    lewis_stress,
    lid_part,
    make_parts,
    planet_states,
    pose,
    print_pose,
    thin_walls,
    wall_thicknesses,
)

PRINTED = {"sun_adapter", "sun", "planet_0", "planet_1", "planet_2", "carrier_output",
           "housing", "base", "lid"}  # fmt: skip


def gear_gap(gear, planet, c: complex, r_gear_side: float, r_planet_side: float, ring: bool):
    """BREP distance gear (sun or housing ring) to a planet centred at c, faces pre-filtered."""
    a = _faces(gear, c, r_gear_side + MARGIN)
    b = _faces(planet, 0j, r_planet_side + (-MARGIN if ring else MARGIN), outside=ring)
    d = BRepExtrema_DistShapeShape()
    d.LoadS1(a.wrapped)
    d.LoadS2(b.wrapped)
    d.SetMultiThread(True)
    d.Perform()
    assert d.IsDone()
    return d.Value()


@pytest.fixture(scope="module")
def home():
    parts = make_parts(Params())
    ctypes.CDLL("libc.so.6").malloc_trim(0)
    return parts


def test_ratio_and_assembly_condition():
    p = Params()
    assert p.ratio == pytest.approx(10.0)
    assert (p.zs * p.z2 + p.zr * p.z1) % p.n_planets == 0
    assert any("assembly" in w for w in design_warnings(Params(zr=55)))


def test_no_warnings_and_walls():
    p = Params()
    assert design_warnings(p) == []
    assert thin_walls(p) == {}
    assert all(v >= p.min_wall - 1e-9 for v in wall_thicknesses(p).values())
    assert "planet small root to 685 seat" in thin_walls(Params(x2=0.2))


def test_lewis_all_gears():
    p = Params()
    s = lewis_stress(p, p.torque_peak)
    assert set(s) == {"sun", "planet_large", "planet_small", "ring"}
    assert max(s.values()) <= p.sigma_allow
    # ring-side force from the planet moment balance: Ft2 / Ft1 = r_large / r_small
    ratio = (p.m1 * p.z1) / (p.m2 * p.z2)
    y14 = 0.154 - 0.912 / 14
    assert s["planet_small"] * p.b2 * y14 * math.pi * p.m2 == pytest.approx(
        ratio * s["sun"] * p.b1 * (0.154 - 0.912 / 15) * math.pi * p.m1
    )


def test_outer_diameter(home):
    p = Params()
    assert 2 * p.R_body <= 84.0
    for n in ("housing", "base", "lid"):
        bb = home[n].bounding_box()
        assert max(bb.size.X, bb.size.Y) <= 84.0, n


def test_bearing_and_pin_fit(home):
    p, z = Params(), Params().zs_
    d, D, B = BEARINGS[p.bearing]
    assert (d, D, B) == (5, 11, 5)
    planet = home["planet_0"]
    c = planet_states(p, 0.0)[0][0]
    # 685 seats (D + 2 fit_clr) at both ends, pin bore of the bearings is the 5 mm pin
    for zz in (z["s0"] + B / 2, z["L1"] - B / 2):
        assert not planet.is_inside((c.real + D / 2 + p.fit_clr / 2, c.imag, zz))
        assert planet.is_inside((c.real + D / 2 + p.fit_clr + 0.3, c.imag, zz))
    assert p.pin_socket >= 8.0
    carrier = home["carrier_output"]
    zc = z["L1"] + p.pin_socket - 0.5
    assert not carrier.is_inside((c.real + p.pin_socket_d / 2 - 0.1, c.imag, zc))
    assert carrier.is_inside((c.real + p.pin_socket_d / 2 + 0.3, c.imag, zc))
    # shoulder bolt: Ø5 shoulder s0 .. socket bottom, head below, M4 end in the tapped hole
    bb = home["_pin_shoulder_bolt_5x30_0"].bounding_box()
    assert bb.size.X == pytest.approx(p.pin_head_d)
    assert bb.min.Z == pytest.approx(z["s0"] - p.pin_head_h)
    assert bb.max.Z == pytest.approx(z["L1"] + p.pin_socket + p.pin_thread_len)
    assert z["s0"] + p.pin_len == pytest.approx(z["L1"] + p.pin_socket)  # bearings on shoulder
    zt = z["L1"] + p.pin_socket + 1
    assert not carrier.is_inside((c.real + p.pin_tap_d / 2 - 0.1, c.imag, zt))
    assert carrier.is_inside((c.real + p.pin_tap_d / 2 + 0.3, c.imag, zt))
    # dowel alternative: same socket, plain Ø5 x pin_len
    dow = Params(pin_kind="dowel")
    assert design_warnings(dow) == []


@pytest.mark.parametrize("k,psi", [(0, 0.0), (1, 1.3), (2, 3.1)])
def test_no_interference_and_tooth_gap(home, k, psi):
    p = Params()
    parts = pose(p, home, psi)
    fixed = {n for n, v in parts.items() if v is home[n]}
    assert interference(parts, skip=fixed) == [], psi
    c = planet_states(p, psi)[k][0]
    planet = parts[f"planet_{k}"]
    g_sun = gear_gap(parts["sun"], planet, c, p.r_tip1, p.r_tip_sun, ring=False)
    assert 0 < g_sun < GAP_MAX, g_sun
    g_ring = gear_gap(parts["housing"], planet, c, p.r_tip2, p.r_ring_tip, ring=True)
    assert 0 < g_ring < GAP_MAX, g_ring


def test_fixed_parts(home):
    fixed = {n: v for n, v in home.items() if pose(Params(), home, 1.0)[n] is v}
    assert {"housing", "base", "lid"} <= set(fixed)
    assert interference(fixed) == []


def test_assembly_paths(home):
    """Cluster lowers into the housing; motor needs the separate base; bayonet lid drops in."""
    p, z = Params(), Params().zs_
    h = home["housing"]
    # nothing of the housing inside the large-gear swing circle above the ring
    swing = Pos(0, 0, z["s1"] + 0.01) * Cylinder(p.a + p.r_tip1, z["lid"] - z["s1"], align=MIN)
    assert (h & swing).volume < 1e-3
    # the adapter passes the housing lip, the motor does not pass the ring
    assert p.adapter_r < p.lip_r and p.motor_d / 2 > p.r_ring_tip
    # lid turned back by lock_deg slides down through the entry slots from above
    lid = lid_part(p, -p.lock_deg)
    for dz in (0.0, 1.5, 3.0, 4.5):
        assert (h & (Pos(0, 0, dz) * lid)).volume < 1e-3, dz
    # locked lid cannot be pulled up: lugs hit the slot ceiling
    assert (h & (Pos(0, 0, 0.5) * home["lid"])).volume > 1.0


def test_print_pose(home):
    assert all(v.is_valid for v in home.values())
    for name in PRINTED:
        bb = print_pose(name, home[name]).bounding_box()
        assert bb.min.Z == pytest.approx(0, abs=1e-6)


def test_build_and_reload(home, tmp_path):
    files = build(Params(), tmp_path)
    assert {f[:-4] for f in files if f.endswith(".stl") and "/" not in f} == PRINTED
    assert len([f for f in files if f.startswith("assembled/")]) == len(home)
    assert set(json.loads(files["parts.json"].read_text())["printed"]) == PRINTED
    for name, f in files.items():
        assert f.stat().st_size > 0
        if name.endswith(".step"):
            assert import_step(str(f)).volume > 0
        elif name.endswith(".stl"):
            assert import_stl(str(f)).area > 0
    assert np.isfinite(import_step(str(files["assembly.step"])).volume)


@pytest.mark.parametrize("name", ["carrier_output", "lid"])
def test_lid_output_cone_prints_without_support(home, name):
    """45 deg hub/lid cone: flipped print pose has no downward face flatter than 45 deg there."""
    z = Params().zs_
    part = home[name]
    top = part.bounding_box().max.Z
    region = [f for f in part.faces()
              if z["c1"] - 1e-6 < f.bounding_box().min.Z < top - 1e-6]  # fmt: skip
    assert any(f.geom_type == GeomType.CONE for f in region)
    for f in region:
        nz = -f.normal_at().Z
        assert nz >= -math.cos(math.radians(45)) - 1e-6, (f.geom_type, f.center())
