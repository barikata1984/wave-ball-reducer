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
    PRINT_FLIP,
    Params,
    build,
    design_warnings,
    lewis_stress,
    lid_part,
    make_parts,
    max_bridge,
    planet_states,
    pose,
    print_pose,
    reference_loads,
    size_stepped,
    thin_walls,
    tooth_deg,
    tooth_r,
    wall_thicknesses,
)

PRINTED = {"sun_adapter", "sun", "planet_0", "planet_1", "planet_2", "carrier_output",
           "housing_upper", "housing_lower", "spigot_ring", "base", "lid"}  # fmt: skip
CLUSTER = ("sun", "planet_0", "planet_1", "planet_2", "_bearing_608ZZ_0", "_bearing_608ZZ_1",
           "_bearing_608ZZ_2", "carrier_output")  # fmt: skip
# widest flat ceiling allowed per part in print pose [mm]; the planets bridge the floor under
# the 608 seat (a ring round the journal bore, outer d = seat d - 2 x cone)
BRIDGE_MAX = {"planet_0": 20.5, "planet_1": 20.5, "planet_2": 20.5}


def vol(a, b) -> float:
    """Common volume of two shapes (0 when disjoint)."""
    c = a & b
    return c.volume if c else 0.0


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
    assert (p.zs, p.z1, p.z2, p.zr, p.m1, p.m2) == (12, 26, 13, 54, 1.1, 1.02)
    assert p.a == pytest.approx(20.9) and p.xr == pytest.approx(0.39, abs=0.01)
    assert (p.zs * p.z2 + p.zr * p.z1) % p.n_planets == 0
    assert any("assembly" in w for w in design_warnings(Params(zr=55)))


def test_no_warnings_and_walls():
    p = Params()
    assert design_warnings(p) == []
    assert thin_walls(p) == {}
    assert all(v >= p.min_wall - 1e-9 for v in wall_thicknesses(p).values())
    assert "planet small root to journal bore" in thin_walls(Params(journal_clr=0.7))
    assert any("608" in w for w in design_warnings(Params(b1=8.5)))


def test_stress_is_reference_only():
    """Lewis and pin loads are reported, never a warning, however large the torque."""
    p = Params(torque_peak=50.0)
    assert design_warnings(p) == []
    s = lewis_stress(p, 0.4)
    assert set(s) == {"sun", "planet_large", "planet_small", "ring"}
    r = reference_loads(p, 0.4)
    assert r["608 load [N]"] + r["journal load [N]"] == pytest.approx(r["Ft1 [N]"] + r["Ft2 [N]"])
    assert r["Ft2 [N]"] / r["Ft1 [N]"] == pytest.approx((p.m1 * p.z1) / (p.m2 * p.z2))
    assert all(v > 0 for v in r.values())


def test_size_matches_model(home):
    """size_stepped's numbers are the model's: body flush with the cap, d 79.96, slack >= 0."""
    p, z = Params(), Params().zs_
    s = size_stepped(p)
    assert p.R_body == p.R_cap and s["D_body"] == s["D_max"]
    assert s["D_max"] == pytest.approx(79.96, abs=0.01)  # max OD as before flush_cap
    q = Params(flush_cap=False)
    assert q.R_body < q.R_cap and 2 * q.R_cap == pytest.approx(s["D_max"])
    assert s["height"] == pytest.approx(sum(s["height terms"].values()))
    assert all(v >= -1e-9 for v in s["slack"].values()), s["slack"]
    for n in ("housing_upper", "housing_lower", "base"):
        bb = home[n].bounding_box()
        assert max(bb.size.X, bb.size.Y) <= s["D_body"] + 0.01, n
    bb = home["lid"].bounding_box()
    assert max(bb.size.X, bb.size.Y) == pytest.approx(s["D_max"], abs=0.01)
    lo = min(v.bounding_box().min.Z for v in home.values())
    hi = max(v.bounding_box().max.Z for v in home.values())
    assert lo == pytest.approx(z["bottom"]) and hi == pytest.approx(z["out"])
    assert hi - lo == pytest.approx(s["height"])


def test_printed_parts_are_one_solid(home):
    for n in PRINTED:
        assert home[n].is_valid and len(home[n].solids()) == 1, n


def test_planet_pin_and_608(home):
    """608 in the large gear, carrier pin through its inner race and the small-gear journal."""
    p, z = Params(), Params().zs_
    d, D, B = BEARINGS[p.bearing]
    assert (d, D, B) == (8, 22, 7)
    planet, carrier = home["planet_0"], home["carrier_output"]
    brg = home["_bearing_608ZZ_0"]
    c = planet_states(p, 0.0)[0][0]
    u = c / abs(c)

    def at(r: float, zz: float, dirn: complex = 1j) -> tuple[float, float, float]:
        q = c + r * u * dirn
        return q.real, q.imag, zz

    # seat: D + 2 fit_clr over B at the large end; bearing flush with the large top face
    zm = z["L1"] - B / 2
    assert not planet.is_inside(at(p.seat_r - 0.02, zm)) and planet.is_inside(
        at(p.seat_r + 0.02, zm)
    )
    assert brg.bounding_box().max.Z == pytest.approx(z["L1"])
    assert brg.bounding_box().min.Z == pytest.approx(z["L1"] - B)
    # the outer race rests on the 45 deg cone; under it min_wall of floor joins the gears
    zf = z["L1"] - B - 0.02
    assert planet.is_inside(at(D / 2 - 0.1, zf - p.bearing_shoulder))
    assert planet.is_inside(at(p.journal_r + 0.5, z["L0"] + p.min_wall - 0.05))
    assert not planet.is_inside(at(p.journal_r + 0.5, z["L0"] + p.min_wall + 0.05))
    # pin: d - 2 pin_fit_clr in the inner race, journal_clr in the small-gear bore
    for zz, r_hole in ((zm, d / 2), ((z["s0"] + z["s1"]) / 2, p.journal_r)):
        assert carrier.is_inside(at(p.pin_r - 0.02, zz))
        assert not carrier.is_inside(at(p.pin_r + 0.02, zz))
        hole_part = brg if r_hole == d / 2 else planet
        assert not hole_part.is_inside(at(r_hole - 0.02, zz))
        assert hole_part.is_inside(at(r_hole + 0.02, zz))
    assert p.journal_r - p.pin_r == pytest.approx(p.journal_clr)
    # the pin reaches the small gear's bottom face and no further
    assert carrier.is_inside(at(0, z["s0"] + p.pin_chamfer + 0.05))
    assert not carrier.is_inside(at(0, z["s0"] - 0.05))
    # collar on the inner-race land, seat_gap above it, clear of the outer race
    zc = z["L1"] + p.seat_gap
    assert carrier.is_inside(at(p.seat_d / 2 - 0.1, zc + 0.05))
    assert not carrier.is_inside(at(p.seat_d / 2 - 0.1, zc - 0.05))
    assert not carrier.is_inside(at(p.seat_d / 2 + 0.1, zc + 0.05))
    assert p.pin_r < p.seat_d / 2 <= p.inner_land_d / 2 < D / 2 - p.bearing_shoulder
    assert vol(carrier, brg) < 1e-3 and vol(planet, brg) < 1e-3
    # pushed down by seat_gap + 0.05 the collar meets the inner race
    assert vol(Pos(0, 0, -(p.seat_gap + 0.05)) * carrier, brg) > 1e-3
    assert any("land" in w for w in design_warnings(Params(seat_d=13.0)))


def test_axial_support(home):
    """Sun bottom face rests on the adapter seat; the hex tip clears the socket ceiling."""
    p, z = Params(), Params().zs_
    sun, adapter = home["sun"], home["sun_adapter"]
    assert sun.bounding_box().min.Z == pytest.approx(z["s0"])
    r = (p.sun_seat_r + p.hex_af / 2 + p.hex_clr + p.hex_chamfer) / 2  # inside the contact ring
    assert adapter.is_inside((r, 0, z["s0"] - 0.05))
    assert not adapter.is_inside((r, 0, z["s0"] + 0.05))
    assert sun.is_inside((r, 0, z["s0"] + 0.05))
    assert not adapter.is_inside((p.sun_seat_r + 0.1, 0, z["s0"] - 0.05))
    tip = z["s0"] + p.hex_engage
    assert adapter.bounding_box().max.Z == pytest.approx(tip)
    assert not sun.is_inside((0, 0, tip + p.gear_gap - 0.05))
    assert sun.is_inside((0, 0, tip + p.gear_gap + 0.05))
    assert any("sun seat" in w for w in design_warnings(Params(sun_seat_r=9.5)))
    # planets over the lip with gear_gap; pins end over it as well
    assert z["s0"] - p.lip_top == pytest.approx(p.gear_gap)


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
    g_ring = gear_gap(parts["housing_upper"], planet, c, p.r_tip2, p.r_ring_tip, ring=True)
    assert 0 < g_ring < GAP_MAX, g_ring


def test_fixed_parts(home):
    fixed = {n: v for n, v in home.items() if pose(Params(), home, 1.0)[n] is v}
    assert {"housing_upper", "housing_lower", "base", "lid"} <= set(fixed)
    assert interference(fixed) == []


def test_assembly_paths(home):
    """Carrier pins slide into the planets; the cluster lowers into the housing; cap screws on."""
    p, z = Params(), Params().zs_
    h = home["housing_upper"]
    # nothing of the housing inside the large-gear swing circle above the ring
    swing = Pos(0, 0, z["s1"] + 0.01) * Cylinder(p.a + p.r_tip1, z["lid"] - z["s1"], align=MIN)
    assert vol(h, swing) < 1e-3
    # step 4: the carrier comes down onto the meshed planets, pins straight through
    rest = {n: home[n] for n in CLUSTER if n != "carrier_output"}
    for dz in (2.0, 8.0, 14.0):
        moved = Pos(0, 0, dz) * home["carrier_output"]
        assert all(vol(moved, v) < 1e-3 for v in rest.values()), dz
    # step 5: the whole cluster lowers into the housing (spur small gears slide axially)
    for dz in (1.0, 4.0, 9.0):
        for n in CLUSTER:
            assert vol(Pos(0, 0, dz) * home[n], h) < 1e-3, (n, dz)
    # the adapter passes the housing lip, the motor does not pass the ring
    assert p.adapter_r < p.lip_r and p.motor_d / 2 > p.r_ring_tip
    # the cap screws on along its thread: turned back and lifted by the lead it clears the
    # housing once the ratchet tooth is above the notches; pulled straight up the thread holds
    for turns in (0.35, 1.0, 2.0, 3.5):
        lid = lid_part(p, 360 * turns, turns * p.thread_pitch)
        assert vol(h, lid) < 1e-3, turns
    assert vol(h, (Pos(0, 0, 0.5) * home["lid"])) > 1.0


def test_lid_seat_and_thread(home):
    """Cap stops on the housing top, so the rib-to-carrier gap is set by geometry, not torque."""
    p, z = Params(), Params().zs_
    h, lid, carrier = home["housing_upper"], home["lid"], home["carrier_output"]
    assert h.bounding_box().max.Z == pytest.approx(z["lid"], abs=1e-3)
    skirt_h = (p.thread_turns + 0.5) * p.thread_pitch
    assert lid.bounding_box().min.Z == pytest.approx(z["lid"] - skirt_h, abs=1e-3)
    r = (p.R_shell + p.R_thread_root) / 2  # wall top, away from the arm, between two notches
    for deg in (tooth_deg(p) + 360 * k / p.n_notches + 4.0 for k in (14, 20, 26)):
        x, y = r * math.cos(math.radians(deg)), r * math.sin(math.radians(deg))
        assert h.is_inside((x, y, z["lid"] - 0.05)) and not h.is_inside((x, y, z["lid"] + 0.05))
        assert lid.is_inside((x, y, z["lid"] + 0.05)) and not lid.is_inside((x, y, z["lid"] - 0.05))
    # rib above the carrier plate: gap lid_axial_clr
    zr = z["lid"] - p.lid_rib_h
    assert zr - z["c1"] == pytest.approx(p.lid_axial_clr)
    rr = p.carrier_r + 0.5
    assert lid.is_inside((rr, 0, zr + 0.02)) and not lid.is_inside((rr, 0, zr - 0.02))
    assert carrier.bounding_box().max.Z > zr  # hub passes the cap bore
    # thread: external on the housing, internal in the skirt, radial thread_clr at the crests
    for deg in (0.0, 120.0, 250.0):
        cs, sn = math.cos(math.radians(deg)), math.sin(math.radians(deg))
        for zz in np.linspace(z["lid"] - 6.0, z["lid"] - 1.0, 11):
            r_in = p.R_thread_root + p.thread_clr - 0.03  # inside the cap's ridge crests
            assert not lid.is_inside((r_in * cs, r_in * sn, zz))
            r_out = p.R_crest + 0.03  # outside the housing's ridge crests
            assert not h.is_inside((r_out * cs, r_out * sn, zz))
            assert lid.is_inside(((p.R_cap - 0.1) * cs, (p.R_cap - 0.1) * sn, zz))
            r_h = p.R_thread_root - 0.05  # min_wall behind the thread root
            assert h.is_inside((r_h * cs, r_h * sn, zz))
    # the thread holds the cap only once it is seated: axial play both ways, no contact
    assert interference({"h": h, "lid": lid}) == []
    # pushed down it hits the top face at once; lifted, the flanks meet after thread_clr
    assert vol(h, (Pos(0, 0, -0.05) * lid)) > 0
    assert vol(h, (Pos(0, 0, p.thread_clr - 0.05) * lid)) < 1e-3
    assert vol(h, (Pos(0, 0, p.thread_clr + 0.1) * lid)) > 0
    assert p.R_thread_root - p.R_shell >= p.min_wall - 1e-9
    assert z["thread0"] > z["s1"]


def test_lid_ratchet(home):
    """Tooth sits in a notch when seated; turning the cap back (loosening) hits a notch wall."""
    p, z = Params(), Params().zs_
    h = home["housing_upper"]
    deg = tooth_deg(p)
    r = sum(tooth_r(p)) / 2
    x, y = r * math.cos(math.radians(deg - 1)), r * math.sin(math.radians(deg - 1))
    assert home["lid"].is_inside((x, y, z["lid"] - p.tooth_h / 2))  # tooth below the plate
    assert not h.is_inside((x, y, z["lid"] - p.tooth_h / 2))  # ... inside a notch
    loosen = math.degrees((p.notch_play + 0.3) / r)
    lift = loosen / 360 * p.thread_pitch
    assert vol(h, lid_part(p, loosen, lift)) > 0.01
    # the arm is free: slot inside it, no skirt under it
    for t in (0.3, 0.7):
        a = math.radians(p.arm_at_deg + t * p.arm_deg)
        rr = p.R_cap - p.arm_w - p.arm_slot / 2
        assert not home["lid"].is_inside((rr * math.cos(a), rr * math.sin(a), z["lid"] + 1.0))
        rs = p.R_cap - p.min_wall / 2
        assert not home["lid"].is_inside((rs * math.cos(a), rs * math.sin(a), z["lid"] - 2.0))
    # thread play (thread_clr, 45 deg) turns through whole notch pitches
    play_deg = p.thread_clr / p.thread_pitch * 360
    assert play_deg / (360 / p.n_notches) == pytest.approx(round(play_deg / (360 / p.n_notches)))


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


def test_housing_print_orientation(home):
    """housing_lower prints upside down on its lip top (split face), housing_upper upright on
    the split face: both bed faces are the full split annulus."""
    p = Params()
    for name, flip in (("housing_lower", True), ("housing_upper", False)):
        part = print_pose(name, home[name])
        assert (name in PRINT_FLIP) == flip
        bed = [f for f in part.faces() if f.bounding_box().max.Z < 1e-6]
        assert sum(f.area for f in bed) > 0.5 * math.pi * (p.R_body**2 - p.R_motor_bore**2), name
        assert _steep_area(part) < 0.1 and max_bridge(part) <= 8.5, name


@pytest.mark.parametrize("name", ["carrier_output", "lid"])
def test_lid_output_cone_prints_without_support(home, name):
    """45 deg hub/lid cone: flipped print pose has no downward face flatter than 45 deg there."""
    z = Params().zs_
    part = home[name]
    top = part.bounding_box().max.Z
    r_max = Params().R_shell  # the cone and bore, not the rim (grip grooves, ratchet arm)
    region = [f for f in part.faces()
              if z["c1"] - 1e-6 < f.bounding_box().min.Z < top - 1e-6
              and math.hypot(f.center().X, f.center().Y) < r_max]  # fmt: skip
    assert any(f.geom_type == GeomType.CONE for f in region)
    for f in region:
        nz = -f.normal_at().Z
        assert nz >= -math.cos(math.radians(45)) - 1e-6, (f.geom_type, f.center())


def _steep_area(part):
    """Area of curved faces facing down flatter than 45 deg, z = 0 excluded."""
    steep = 0.0
    us = np.linspace(0.05, 0.95, 10)
    for f in part.faces():
        if f.bounding_box().max.Z < 1e-6 or f.geom_type == GeomType.PLANE:
            continue
        nz = np.array([f.normal_at(u, v).Z for u in us for v in us])
        steep += f.area * np.mean(nz < -math.cos(math.radians(45)) - 1e-3)
    return steep


@pytest.mark.parametrize("name", sorted(PRINTED))
def test_print_without_support(home, name):
    """Print pose: no face flatter than 45 deg facing down; flat ceilings only as bridges no
    wider than 8.5 mm (socket and hole ends), the planets' 608 floor up to BRIDGE_MAX."""
    part = print_pose(name, home[name])
    assert _steep_area(part) < 0.1
    assert max_bridge(part) <= BRIDGE_MAX.get(name, 8.5), max_bridge(part)


def test_housing_split_fit(home):
    """Split at the ring bottom: the lip is in the lower; the loose spigot ring centres both."""
    p, z = Params(), Params().zs_
    up, lo, ring = home["housing_upper"], home["housing_lower"], home["spigot_ring"]
    assert z["split"] == pytest.approx(z["ring_lo"]) and z["split"] == pytest.approx(p.lip_top)
    assert up.bounding_box().min.Z == pytest.approx(z["split"])
    assert lo.bounding_box().max.Z == pytest.approx(z["split"])
    # lip (floor over the motor, bore lip_r) belongs to the lower, the upper starts at the teeth
    zl = (p.motor_top_clr + p.lip_top) / 2
    assert lo.is_inside((p.lip_r + 0.5, 0, zl)) and not lo.is_inside((p.lip_r - 0.5, 0, zl))
    assert not up.is_inside((p.lip_r + 0.5, 0, zl))
    assert interference({"up": up, "lo": lo, "ring": ring}) == []
    r0, r1 = p.spigot_r
    for part, zz in ((up, z["split"] + p.spigot_h / 2), (lo, z["split"] - p.spigot_h / 2)):
        assert ring.is_inside((r1 - 0.02, 0, zz)) and not ring.is_inside((r1 + 0.02, 0, zz))
        assert not part.is_inside((r1 + p.slide_clr - 0.02, 0, zz))
        assert part.is_inside((r1 + p.slide_clr + 0.02, 0, zz))
        assert not part.is_inside((r0 - p.spigot_in_clr + 0.02, 0, zz))
        assert part.is_inside((r0 - p.spigot_in_clr - 0.02, 0, zz))
        assert vol(part, Pos(p.slide_clr + 0.05, 0, 0) * ring) > 0  # radial play limited
    # upper rests on the lower's flat top (pushed down it collides)
    assert vol(lo, (Pos(0, 0, -0.05) * up)) > 0


def test_housing_joint_screws(home):
    """Joint screw holes line up, leave min_wall, and a key reaches them from above."""
    p, z = Params(), Params().zs_
    up, lo = home["housing_upper"], home["housing_lower"]
    L = p.screw_len(p.joint_clamp)
    assert L == 12 and L - p.joint_clamp >= 6  # stock cap screw, 6+ mm in the tube
    # head (3 mm) inside the counterbore below the ring top; hex hole longer than the screw tip
    assert z["joint_seat"] + p.joint_head_h <= z["s1"]
    assert L - p.joint_clamp + p.hex_hole_extra <= z["split"] - z["base_top"]
    for k in range(p.n_joint):
        a = math.radians(p.joint_phase_deg) + 2 * math.pi * k / p.n_joint
        x, y = p.joint_r * math.cos(a), p.joint_r * math.sin(a)
        c, s = math.cos(a), math.sin(a)
        assert not up.is_inside((x, y, z["split"] + 1))
        assert up.is_inside((x + 2.2 * c, y + 2.2 * s, z["split"] + 1))
        r_seat = (p.m3_clear_d / 2 + p.cbore_r) / 2
        seat = (x + r_seat * c, y + r_seat * s, z["joint_seat"])
        assert up.is_inside((*seat[:2], seat[2] - 0.05))
        assert not up.is_inside((*seat[:2], seat[2] + 0.05))
        assert not lo.is_inside((x, y, z["split"] - (L - p.joint_clamp) + 0.1))
        screw = Pos(x, y, z["joint_seat"] - L) * Cylinder(p.screw_d / 2 - 0.3, L, align=MIN)
        assert vol(screw, up) < 1e-3 and vol(screw, lo) < 1e-3
        col = Pos(x, y, z["joint_seat"] + 0.01) * Cylinder(p.cbore_r - 0.02, z["lid"], align=MIN)
        assert vol(col, up) < 1e-3
        assert p.joint_r + p.cbore_r <= p.R_shell + 1e-9
    walls = wall_thicknesses(p)
    for k in ("upper: joint counterbore to ring root", "upper: joint counterbore to outside",
              "lower: joint hole to motor bore", "lower: joint hole to outside",
              "spigot groove to joint hole", "spigot groove to outside"):  # fmt: skip
        assert walls[k] >= p.min_wall - 1e-9, k


def test_assembly_order(home):
    """Lower over motor + adapter, upper over the adapter, screws in, then cluster, then cap."""
    p = Params()
    up, lo = home["housing_upper"], home["housing_lower"]
    assert p.R_motor_bore > p.motor_d / 2 and p.R_motor_bore > p.adapter_r
    assert p.lip_r > p.adapter_r and up.bounding_box().min.Z >= 0
    # the lip in the lower passes the adapter: lower goes on over motor + adapter in one move
    assert vol(lo, home["sun_adapter"]) < 1e-3
    for n in (n for n in home if n != "housing_upper"):
        if n in CLUSTER or n == "lid":
            continue  # cluster and cap go in after the upper is screwed down
        assert vol(up, home[n]) < 1e-3, n
    assert vol(lo, home["_motor_GM5208-24"]) < 1e-3
    for dz in (30.0, 10.0, p.spigot_h + 0.1):
        assert vol(lo, (Pos(0, 0, dz) * up)) < 1e-3, dz
        assert vol(home["spigot_ring"], (Pos(0, 0, dz) * up)) < 1e-3, dz


def test_m8_bolt_variant():
    """pin_kind='m8_bolt': bolts through the carrier hub into the planets, nothing collides."""
    p = Params(pin_kind="m8_bolt", n_out_holes=3, out_phase_deg=60.0)  # holes between heads
    assert any("output hole" in w for w in design_warnings(Params(pin_kind="m8_bolt")))
    assert design_warnings(p) == [] and thin_walls(p) == {}
    parts = make_parts(p)
    bolts = [n for n in parts if n.startswith("_bolt_M8")]
    assert len(bolts) == p.n_planets and p.bolt_len in (16, 20, 25)
    for psi in (0.0, 2.0):
        q = pose(p, parts, psi)
        sub = {n: q[n] for n in (*bolts, "planet_0", "planet_1", "planet_2", "_bearing_608ZZ_0",
                                 "lid", "housing_upper", "sun")}  # fmt: skip
        assert interference(sub) == [], psi
    # the bolt cuts its thread in the carrier: overlap only between tap and nominal diameter
    d, L = 8.0, p.bolt_len
    cut = vol(parts[bolts[0]], parts["carrier_output"])
    assert 0 < cut <= math.pi / 4 * (d**2 - p.bolt_tap_d**2) * L + 1e-6
