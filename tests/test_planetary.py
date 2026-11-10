import json
import math

import numpy as np
import pytest
from build123d import import_step, import_stl

from planetary import (
    BUILD_VOLUME,
    Params,
    build,
    design_warnings,
    lewis_stress,
    make_parts,
    planet_states,
    pose,
    print_pose,
    thin_walls,
)
from reducer import interference

PRINTED = {"sun_adapter", "sun", "planet_0", "planet_1", "planet_2", "carrier_output",
           "housing", "ring_lid"}  # fmt: skip
# BREP distance on helical solids takes ~45 s per pair: one planet per pose, a different one each
GEAR_PAIRS = [("sun", "planet_{}"), ("housing", "planet_{}"), ("ring_lid", "planet_{}")]


@pytest.fixture(scope="module")
def home():
    return make_parts(Params())


def test_tooth_conditions():
    p = Params()
    assert p.zr == p.zs + 2 * p.zp
    assert (p.zs + p.zr) % p.n_planets == 0
    assert p.zp + 2 < (p.zs + p.zp) * math.sin(math.pi / p.n_planets)
    assert design_warnings(p) == [] or all("tip land" in w for w in design_warnings(p))
    assert any("concentric" in w for w in design_warnings(Params(zr=61)))
    assert any("assembly" in w for w in design_warnings(Params(zs=13, zr=61, zp=24)))
    assert any("undercut" in w for w in design_warnings(Params(xs=0.0, xp=0.0)))
    assert any("KHK" in w for w in design_warnings(Params(zs=36, zp=12, zr=60)))


def test_ratio():
    p = Params()
    assert p.ratio == 1 + p.zr / p.zs == 6


def test_carrier_speed_from_rolling():
    """Planet pitch points roll without slip on the fixed ring and on the sun."""
    p = Params()
    rp, rs, rr = p.r(p.zp), p.r(p.zs), p.r(p.zr)
    d = 1e-6
    for psi in (0.0, 1.0, 4.0):
        for (c0, t0), (c1, t1) in zip(planet_states(p, psi), planet_states(p, psi + d)):
            w = (t1 - t0) / d  # planet spin per unit input
            v = (c1 - c0) / d  # planet centre velocity
            u = c0 / abs(c0)  # outward unit vector
            at_ring = v + 1j * w * rp * u  # velocity of the planet point touching the ring
            at_sun = v - 1j * w * rp * u
            assert abs(at_ring) < 1e-6
            assert at_sun == pytest.approx(1j * rs * u, abs=1e-6)  # sun surface speed
            assert abs(v) == pytest.approx((rr - rp) / p.ratio, rel=1e-6)


def test_lewis():
    p = Params()
    s = lewis_stress(p, 0.25)
    assert s["sun"] > s["planet"] > 0
    assert lewis_stress(p, 0.5)["sun"] == pytest.approx(2 * s["sun"])
    assert max(lewis_stress(p, p.torque_peak).values()) < p.sigma_allow


def test_min_wall():
    assert thin_walls(Params()) == {}
    assert "planet root to bearing seat" in thin_walls(Params(m=1.1))


@pytest.mark.parametrize("k,psi", [(0, 0.0), (1, 0.9), (2, 2.3)])
def test_no_interference_and_tooth_gap(home, k, psi):
    p = Params()
    parts = pose(p, home, psi)
    assert interference(parts) == [], psi
    for a, b in GEAR_PAIRS:
        gap = parts[a].distance_to(parts[b.format(k)])
        assert 0 < gap < 0.3, (a, b, gap)  # meshing, flank clearance only


def test_print_pose_and_build_volume(home):
    for name in PRINTED:
        bb = print_pose(name, home[name]).bounding_box()
        assert bb.min.Z == pytest.approx(0, abs=1e-6)
        assert bb.size.X < BUILD_VOLUME[0] and bb.size.Y < BUILD_VOLUME[1]
        assert bb.size.Z < BUILD_VOLUME[2]


def test_build_and_reload(tmp_path):
    files = build(Params(), tmp_path)
    assert {f[:-4] for f in files if f.endswith(".stl") and "/" not in f} == PRINTED
    assert len([f for f in files if f.startswith("assembled/")]) == len(make_parts(Params()))
    assert set(json.loads(files["parts.json"].read_text())["printed"]) == PRINTED
    for name, f in files.items():
        assert f.stat().st_size > 0
        if name.endswith(".step"):
            assert import_step(str(f)).volume > 0
        elif name.endswith(".stl"):
            assert import_stl(str(f)).area > 0
    assert np.isfinite(import_step(str(files["assembly.step"])).volume)
