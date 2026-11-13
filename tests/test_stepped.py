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
    tooth_deg,
    wall_thicknesses,
)

PRINTED = {"sun_adapter", "sun", "planet_0", "planet_1", "planet_2", "carrier_output",
           "housing_upper", "housing_lower", "base", "lid", "pin_link"}  # fmt: skip


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
    for n in ("housing_upper", "housing_lower", "base", "lid"):
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
    assert p.pin_socket >= p.pin_socket_min - 1e-9
    carrier = home["carrier_output"]
    zc = z["L1"] + p.pin_socket - 0.5
    assert not carrier.is_inside((c.real + p.pin_socket_d / 2 - 0.1, c.imag, zc))
    assert carrier.is_inside((c.real + p.pin_socket_d / 2 + 0.3, c.imag, zc))
    # shoulder bolt: Ø5 shoulder s0 .. socket bottom, head below, M4 end in the tapped hole
    bb = home["_pin_shoulder_bolt_5x30_0"].bounding_box()
    assert bb.size.X == pytest.approx(p.pin_head_d)
    assert bb.min.Z == pytest.approx(p.head_top - p.pin_head_h)
    assert bb.max.Z == pytest.approx(z["L1"] + p.pin_socket + p.pin_thread_len)
    assert p.head_top + p.pin_len == pytest.approx(z["L1"] + p.pin_socket)
    zt = z["L1"] + p.pin_socket + 1
    assert not carrier.is_inside((c.real + p.pin_tap_d / 2 - 0.1, c.imag, zt))
    assert carrier.is_inside((c.real + p.pin_tap_d / 2 + 0.3, c.imag, zt))
    # dowel alternative: same socket, plain Ø5 x pin_len
    dow = Params(pin_kind="dowel")
    assert design_warnings(dow) == []


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
    c = planet_states(p, 0.0)[0][0]
    head = home["_pin_shoulder_bolt_5x30_0"]
    assert head.is_inside((c.real + p.pin_head_d / 2 - 0.1, c.imag, p.head_top - 0.05))
    assert not head.is_inside((c.real + 3.0, c.imag, p.head_top + 0.05))
    carrier, zs = home["carrier_output"], z["L1"] + p.seat_gap
    assert p.seat_d <= 6.2 and p.seat_gap > 0 and p.pin_head_gap > 0
    assert carrier.is_inside((c.real + p.seat_d / 2 - 0.1, c.imag, zs + 0.05))
    assert not carrier.is_inside((c.real + p.seat_d / 2 - 0.1, c.imag, zs - 0.05))
    assert not carrier.is_inside((c.real + p.seat_d / 2 + 0.1, c.imag, zs + 0.05))


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
    """Cluster lowers into the housing; motor needs the separate base; the lid screws in."""
    p, z = Params(), Params().zs_
    h = home["housing_upper"]
    # nothing of the housing inside the large-gear swing circle above the ring
    swing = Pos(0, 0, z["s1"] + 0.01) * Cylinder(p.a + p.r_tip1, z["lid"] - z["s1"], align=MIN)
    assert (h & swing).volume < 1e-3
    # the adapter passes the housing lip, the motor does not pass the ring
    assert p.adapter_r < p.lip_r and p.motor_d / 2 > p.r_ring_tip
    # the lid screws in along its thread: turned back and lifted by the lead it clears the
    # housing once the ratchet tooth is above the notches; pulled straight up the thread holds
    for turns in (0.35, 1.0, 2.0, 3.5):
        lid = lid_part(p, 360 * turns, turns * p.thread_pitch)
        assert (h & lid).volume < 1e-3, turns
    assert (h & (Pos(0, 0, 0.5) * home["lid"])).volume > 1.0


def test_lid_seat_and_thread(home):
    """Lid stops on the housing top, so the rib-to-carrier gap is set by geometry, not by torque."""
    p, z = Params(), Params().zs_
    h, lid, carrier = home["housing_upper"], home["lid"], home["carrier_output"]
    assert h.bounding_box().max.Z == pytest.approx(z["lid"])
    assert lid.bounding_box().min.Z == pytest.approx(z["band"])
    r = p.R_body - p.grip_d - 0.3  # shoulder ring, away from the arm, between two notches
    for deg in (tooth_deg(p) + 360 * k / p.n_notches + 4.0 for k in (14, 20, 26)):
        x, y = r * math.cos(math.radians(deg)), r * math.sin(math.radians(deg))
        assert h.is_inside((x, y, z["lid"] - 0.05)) and not h.is_inside((x, y, z["lid"] + 0.05))
        assert lid.is_inside((x, y, z["lid"] + 0.05)) and not lid.is_inside((x, y, z["lid"] - 0.05))
    # rib above the carrier plate: gap lid_axial_clr
    zr = z["lid"] - p.lid_rib_h
    assert zr - z["c1"] == pytest.approx(p.lid_axial_clr)
    rr = p.carrier_r + 0.5
    assert lid.is_inside((rr, 0, zr + 0.02)) and not lid.is_inside((rr, 0, zr - 0.02))
    assert carrier.bounding_box().max.Z > zr  # hub passes the lid bore
    # the thread holds the lid only once it is seated: axial play both ways, no contact
    assert interference({"h": h, "lid": lid}) == []
    # pushed down it hits the shoulder at once; lifted, the flanks meet after thread_clr (45 deg)
    assert (h & (Pos(0, 0, -0.05) * lid)).volume > 0
    assert (h & (Pos(0, 0, p.thread_clr - 0.05) * lid)).volume < 1e-3
    assert (h & (Pos(0, 0, p.thread_clr + 0.1) * lid)).volume > 0
    # shell behind the groove root keeps min_wall over the whole thread zone
    rg = p.R_sweep + p.thread_depth
    zn = z["lid"] - p.tooth_h - 0.3  # below the ratchet notches
    for zz in np.linspace(z["lid"] - (p.thread_turns + 1) * p.thread_pitch, zn, 25):
        for deg in (0.0, 90.0, 180.0):
            c, s = math.cos(math.radians(deg)), math.sin(math.radians(deg))
            assert h.is_inside(((rg + 0.02) * c, (rg + 0.02) * s, zz)), (zz, deg)
    assert p.R_body - rg >= p.min_wall - 1e-9
    assert any("groove" in w for w in design_warnings(Params(thread_turns=6)))


def test_lid_ratchet(home):
    """Tooth sits in a notch when seated; turning the lid back (loosening) runs into a notch wall."""
    p, z = Params(), Params().zs_
    h = home["housing_upper"]
    deg = tooth_deg(p)
    r = (p.R_sweep + p.R_body) / 2
    x, y = r * math.cos(math.radians(deg - 1)), r * math.sin(math.radians(deg - 1))
    assert home["lid"].is_inside((x, y, z["lid"] - p.tooth_h / 2))  # tooth below the lid face
    assert not h.is_inside((x, y, z["lid"] - p.tooth_h / 2))  # ... inside a notch
    loosen = math.degrees((p.notch_play + 0.3) / r)
    lift = loosen / 360 * p.thread_pitch
    assert (h & lid_part(p, loosen, lift)).volume > 0.01
    # the arm is free: nothing of the lid between arm and band over the arm span
    for t in (0.3, 0.7):
        a = math.radians(p.arm_at_deg + t * p.arm_deg)
        rr = p.R_body - p.arm_w - p.arm_slot / 2
        assert not home["lid"].is_inside((rr * math.cos(a), rr * math.sin(a), z["lid"] + 1.0))
        assert not home["lid"].is_inside(
            ((rr + 0.5) * math.cos(a), (rr + 0.5) * math.sin(a), z["lid"] - 0.3)
        )
    # 40 notches: the thread's axial play (thread_clr, 45 deg) turns through whole notch pitches
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


def _down_faces(part):
    """(planar ceilings, area of curved faces facing down flatter than 45 deg), z = 0 excluded."""
    flats, steep = [], 0.0
    us = np.linspace(0.05, 0.95, 10)
    for f in part.faces():
        if f.bounding_box().max.Z < 1e-6:
            continue
        if f.geom_type == GeomType.PLANE:
            if f.normal_at().Z < -0.999 and f.area > 0.1:
                flats.append(f)
            continue
        nz = np.array([f.normal_at(u, v).Z for u in us for v in us])
        steep += f.area * np.mean(nz < -math.cos(math.radians(45)) - 1e-3)
    return flats, steep


@pytest.mark.parametrize("name", sorted(PRINTED))
def test_print_without_support(home, name):
    """Print pose: no face flatter than 45 deg facing down, ceilings only short bridges (<= 8 mm
    across: socket ends, tapped-hole ends). No exception: the split housing needs no support."""
    part = print_pose(name, home[name])
    flats, steep = _down_faces(part)
    assert steep < 0.1, steep
    for f in flats:
        bb = f.bounding_box()
        assert min(bb.size.X, bb.size.Y) <= 8.5, (name, bb.min.Z, f.area)


def test_housing_split_fit(home):
    """Spigot centres the upper on its outer wall (radial slide_clr); inner wall and roof clear."""
    p, z = Params(), Params().zs_
    up, lo = home["housing_upper"], home["housing_lower"]
    assert up.bounding_box().min.Z == pytest.approx(z["split"])
    assert lo.bounding_box().max.Z == pytest.approx(z["split"] + p.spigot_h)
    assert interference({"up": up, "lo": lo}) == []
    r0, r1 = p.spigot_r
    zz = z["split"] + p.spigot_h / 2
    assert lo.is_inside((r1 - 0.02, 0, zz)) and not lo.is_inside((r1 + 0.02, 0, zz))
    assert not up.is_inside((r1 + p.slide_clr - 0.02, 0, zz))
    assert up.is_inside((r1 + p.slide_clr + 0.02, 0, zz))
    assert not up.is_inside((r0 - p.spigot_in_clr + 0.02, 0, zz))
    assert up.is_inside((r0 - p.spigot_in_clr - 0.02, 0, zz))
    # centring: worst-case offset of the ring = radial clearance of the locating wall
    assert p.slide_clr <= 0.05
    # upper rests on the tube top (pushed down it collides), spigot engaged before contact
    assert (lo & (Pos(0, 0, -0.05) * up)).volume > 0
    assert (up & (Pos(p.slide_clr + 0.05, 0, 0) * lo)).volume > 0  # radial play limited
    assert p.r_ring_tip + p.min_wall <= r0 - p.spigot_in_clr + 1e-9
    assert any("overhangs" in w for w in design_warnings(Params(spigot_in_clr=-0.1)))


def test_housing_joint_screws(home):
    """Joint screw holes line up, leave min_wall, and a key reaches them from above."""
    p, z = Params(), Params().zs_
    up, lo = home["housing_upper"], home["housing_lower"]
    L = p.screw_len(p.joint_clamp)
    assert L in (6, 8, 12, 16) and L - p.joint_clamp >= 6  # stock cap screw, 6+ mm in the tube
    for k in range(p.n_joint):
        a = math.radians(p.joint_phase_deg) + 2 * math.pi * k / p.n_joint
        x, y = p.joint_r * math.cos(a), p.joint_r * math.sin(a)
        # clearance hole through the slab, head seat at joint_seat, hex hole in the tube
        c, s = math.cos(a), math.sin(a)
        assert not up.is_inside((x, y, z["split"] + 1))
        assert up.is_inside((x + 2.2 * c, y + 2.2 * s, z["split"] + 1))
        r_seat = (p.m3_clear_d / 2 + p.cbore_r) / 2
        seat = (x + r_seat * c, y + r_seat * s, z["joint_seat"])
        assert up.is_inside((*seat[:2], seat[2] - 0.05))
        assert not up.is_inside((*seat[:2], seat[2] + 0.05))
        assert not lo.is_inside((x, y, z["split"] - (L - p.joint_clamp) + 0.1))
        # the screw (shank d 3 down to L below the seat, head cbore above) hits neither part
        screw = Pos(x, y, z["joint_seat"] - L) * Cylinder(p.screw_d / 2 - 0.3, L, align=MIN)
        assert (screw & up).volume < 1e-3 and (screw & lo).volume < 1e-3
        # driver access: nothing of the upper above the head seat over the counterbore
        col = Pos(x, y, z["joint_seat"] + 0.01) * Cylinder(p.cbore_r - 0.02, z["lid"], align=MIN)
        assert (col & up).volume < 1e-3
        # counterbore inside the swing circle is a hole only: no boss there
        assert p.joint_r + p.cbore_r < p.R_sweep
    walls = wall_thicknesses(p)
    for k in ("upper: joint counterbore to ring root", "upper: joint counterbore to outside",
              "lower: joint hole to motor bore", "lower: joint hole to outside",
              "upper: spigot groove to ring bore"):  # fmt: skip
        assert walls[k] >= p.min_wall - 1e-9, k


def test_assembly_order(home):
    """Lower over motor + adapter, upper over the adapter, screws in, then cluster, then lid."""
    p, z = Params(), Params().zs_
    up, lo = home["housing_upper"], home["housing_lower"]
    # the lower slides down over the motor and the adapter
    assert p.R_motor_bore > p.motor_d / 2 and p.R_motor_bore > p.adapter_r
    # the upper's lip clears the turning adapter, and the upper stays above the motor
    assert p.lip_r > p.adapter_r and up.bounding_box().min.Z >= 0
    motor = home["_motor_GM5208-24"]
    for n in (n for n in home if n != "housing_upper"):
        if n.startswith(("planet", "_bearing", "_pin", "pin_link", "carrier", "sun", "lid")):
            continue  # cluster and lid go in after the upper is screwed down
        assert (up & home[n]).volume < 1e-3, n
    assert (lo & motor).volume < 1e-3
    # lowering the upper by its full height from above never meets the lower except at the seat
    for dz in (30.0, 10.0, p.spigot_h + 0.1):
        assert (lo & (Pos(0, 0, dz) * up)).volume < 1e-3, dz
    # ultra-low rotor screw heads: s0 now set by the lip + shoulder-bolt-head stack
    assert z["s0"] == pytest.approx(
        p.motor_top_clr + p.min_wall + p.gear_gap + p.pin_head_h + p.link_t + p.pin_head_gap
    )
    assert z["s0"] >= p.adapter_t + p.head_h + p.head_gap


def test_pin_link(home):
    """Link hangs on the 3 bolt heads, fits the shoulders, clears everything that turns against
    it, and goes in with the cluster through the ring bore."""
    p, z = Params(), Params().zs_
    link = home["pin_link"]
    d = BEARINGS[p.bearing][0]
    r0, r1 = p.link_r
    bb = link.bounding_box()
    assert bb.min.Z == pytest.approx(p.head_top) and bb.size.Z == pytest.approx(p.link_t)
    assert p.link_t >= p.min_wall
    for k, (c, _) in enumerate(planet_states(p, 0.0)):
        head = home[f"_pin_shoulder_bolt_5x30_{k}"]
        u = c / abs(c)
        # hole: shoulder d + 2 slide_clr; link rests on the head top (touching, no overlap)
        for r, inside in ((d / 2 + p.slide_clr - 0.02, False), (d / 2 + p.slide_clr + 0.02, True)):
            q = c + r * u * 1j
            assert link.is_inside((q.real, q.imag, p.head_top + p.link_t / 2)) is inside
        q = c + 3.5 * u * 1j
        assert head.is_inside((q.real, q.imag, p.head_top - 0.05))
        assert link.is_inside((q.real, q.imag, p.head_top + 0.05))
        assert (link & head).volume < 1e-3
    # axial gaps: link top to the planets / small bearings, heads to the housing lip
    zt = p.head_top + p.link_t
    assert z["s0"] - zt == pytest.approx(p.pin_head_gap)
    assert p.head_top - p.pin_head_h - p.lip_top == pytest.approx(p.gear_gap)
    # radial gaps: adapter (+ rotor screw heads) inside, housing ring bore outside
    assert r0 - p.adapter_r >= p.gear_gap and p.r_ring_tip - r1 >= p.gear_gap
    assert p.head_top >= p.adapter_t + p.head_h + p.gear_gap  # above the rotor screw heads
    h = home["housing_upper"]
    for a in np.linspace(0, 2 * math.pi, 7):
        rr = r1 + p.gear_gap
        assert not h.is_inside((rr * math.cos(a), rr * math.sin(a), zt - 0.05))
    # nothing that turns against the link comes within gear_gap of it (any input angle)
    for psi in (0.0, 1.3, 3.1):
        parts = pose(p, home, psi)
        turning = ["sun", "sun_adapter", "housing_upper", "planet_0", "planet_1", "planet_2"]
        turning += [n for n in parts if n.startswith("_bearing") and "_small_" in n]
        for n in turning:
            for dz in (p.gear_gap - 0.05, -(p.gear_gap - 0.05)):
                assert (Pos(0, 0, dz) * parts["pin_link"] & parts[n]).volume < 1e-3, (n, psi, dz)
        # ... and moved by just over gear_gap it does meet the planets' bearings
        moved = Pos(0, 0, p.gear_gap + 0.05) * parts["pin_link"]
        assert any((moved & parts[n]).volume > 1e-3 for n in turning)
    # the cluster (link included) passes the ring bore from above
    assert r1 < p.r_ring_tip
    assert any("pin link" in w for w in design_warnings(Params(link_w=14.0)))
