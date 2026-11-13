"""Unit tests for components.py on small stand-alone geometry (no reducer)."""

import math

import pytest
from build123d import Cylinder, Pos

import components as cp
from components import MIN

R, DEPTH, P, CLR = 20.0, 0.8, 2.0, 0.3  # thread root radius, depth, pitch, radial clearance
ZS = 10.0  # cap seat / housing top


def _band() -> object:
    """Housing top band: tube R-1.2 .. R up to ZS with housing_thread on it."""
    tube = Cylinder(R + cp.ROOT_BURY, ZS, align=MIN) - Cylinder(R - 1.2, ZS, align=MIN)
    return tube + cp.housing_thread(R, DEPTH, P, ZS, 4)


@pytest.fixture(scope="module")
def band():
    return _band()


@pytest.fixture(scope="module")
def cap():
    return cp.screw_cap(R, DEPTH, P, CLR, 1.2, ZS, 2.0, 3)


def test_sector_volume():
    s = cp.sector(5, 7, 10, 100, 1, 2)
    assert s.volume == pytest.approx(math.pi * (49 - 25) * 2 * 90 / 360, rel=1e-6)


def test_thread_45_profile_and_hand():
    t = cp.thread_45(R, DEPTH, P, 0.0, 6.0)
    assert t.is_valid
    # ridge centres at angle 0: z0 + P/2 + k P; at 90 deg a quarter pitch higher (right hand)
    r = R + DEPTH - 0.05
    for zc in (1.0, 3.0, 5.0):
        assert t.is_inside((r, 0, zc)) and not t.is_inside((r, 0, zc + P / 2))
        assert t.is_inside((0, r, zc + P / 4))
    # 45 deg flanks: width at mid-depth = (root + crest)/2 = depth + 0.2
    rm = R + DEPTH / 2
    w = sum(0.01 for k in range(200) if t.is_inside((rm, 0, 2.0 + k * 0.01)))
    assert w == pytest.approx(DEPTH + 0.2, abs=0.03)
    assert t.bounding_box().max.Z == pytest.approx(6.0, abs=1e-3)


def test_cap_seat_and_thread_play(band, cap):
    assert cap.bounding_box().min.Z == pytest.approx(ZS - 3.5 * P, abs=1e-3)
    assert max(cap.bounding_box().size.X, cap.bounding_box().size.Y) == pytest.approx(
        2 * cp.cap_outer_r(R, DEPTH, CLR, 1.2), abs=1e-3
    )
    assert (band & cap).volume < 1e-3  # seated: ridges mid-groove, no contact
    assert (band & (Pos(0, 0, -0.05) * cap)).volume > 0  # plate on the top face
    assert (band & (Pos(0, 0, CLR - 0.05) * cap)).volume < 1e-3  # axial play = CLR
    assert (band & (Pos(0, 0, CLR + 0.1) * cap)).volume > 0  # then the flanks meet
    for turns in (0.25, 1.0, 2.5):  # screwing off along the helix stays free
        moved = Pos(0, 0, turns * P) * cap.rotate(cp.Axis.Z, 360 * turns)
        assert (band & moved).volume < 1e-3, turns


def test_ratchet_arm_and_notches(band, cap):
    deg, rm = 30.0, R - 0.6
    arm = cp.ratchet_arm(cap, cp.cap_outer_r(R, DEPTH, CLR, 1.2), 3.5, 0.8, 20.0, 40.0, ZS, 2.0,
                         3.5 * P, deg, (R - 1.1, R - 0.1), 0.6, 1.0)  # fmt: skip
    notched = band - cp.ratchet_notches_cut(R - 2, R + 2, ZS, 36, deg, rm, 0.6, 1.0, 0.4)
    assert (notched & arm).volume < 1e-3
    x, y = rm * math.cos(math.radians(deg - 1)), rm * math.sin(math.radians(deg - 1))
    assert arm.is_inside((x, y, ZS - 0.3)) and not notched.is_inside((x, y, ZS - 0.3))
    loosen = math.degrees(0.7 / rm)
    moved = Pos(0, 0, loosen / 360 * P) * arm.rotate(cp.Axis.Z, loosen)
    assert (notched & moved).volume > 0.01  # the vertical face hits the notch wall
    a = math.radians(20.0 + 20.0)
    ro = cp.cap_outer_r(R, DEPTH, CLR, 1.2) - 0.3
    assert not arm.is_inside((ro * math.cos(a), ro * math.sin(a), ZS - 2.0))  # no skirt under arm


def test_grip_grooves():
    cuts = cp.grip_grooves_cut(10.0, 0.0, 2.0, 12, 1.0, 0.5)
    disc = Cylinder(10.0, 2.0, align=MIN) - cuts
    assert len(cuts) == 12
    assert disc.volume == pytest.approx(math.pi * 100 * 2 - 12 * 1.0 * 0.5 * 2, rel=2e-3)


def test_split_housing():
    tube = Pos(0, 0, -10) * (Cylinder(32, 25, align=MIN) - Cylinder(22, 25, align=MIN))
    up, lo, ring = cp.split_housing(tube, 0.0, 32.0, 29.0, 1.6, 1.5, 0.3, 0.05, 4, 25.5, 45.0,
                                    4.0, 12.0, 3.0, 3.4, 5.5, 12)  # fmt: skip
    # both split faces flat at z = 0 (bed faces); the ring sits half in each groove
    assert up.bounding_box().min.Z == pytest.approx(0.0)
    assert lo.bounding_box().max.Z == pytest.approx(0.0)
    assert ring.bounding_box().min.Z == pytest.approx(-1.5)
    assert ring.bounding_box().max.Z == pytest.approx(1.5)
    assert (up & lo).volume < 1e-3 and (up & ring).volume < 1e-3 and (lo & ring).volume < 1e-3
    for part, zz in ((up, 0.75), (lo, -0.75)):
        assert not part.is_inside((30.6 + 0.03, 0, zz)) and part.is_inside((30.65 + 0.03, 0, zz))
        assert not part.is_inside((28.7 + 0.02, 0, zz)) and part.is_inside((28.7 - 0.02, 0, zz))
    x, y = 25.5 * math.cos(math.radians(45)), 25.5 * math.sin(math.radians(45))
    assert not up.is_inside((x, y, 2.0)) and not up.is_inside((x + 1.5, y, 6.0))  # hole, cbore
    assert not lo.is_inside((x, y, -8.9)) and lo.is_inside((x, y, -9.1))  # hex hole 9 deep


def test_sun_adapter_and_socket():
    a = cp.sun_adapter(
        16.0, 3.0, 8.0, 5.0, 7.0, 0.15, 4.5, 0.4, 4, 25.0, 3.4, (6.2, 4.0), -20.0, -30.0, 3.2, 2.1
    )
    assert a.bounding_box().max.Z == pytest.approx(5.0)
    assert a.bounding_box().min.Z == pytest.approx(-30.0)
    assert a.is_inside((7.5, 0, 4.9)) and not a.is_inside((7.5, 0, 5.1))  # seat top
    assert not a.is_inside((0, 3.6, 2.0)) and a.is_inside((0, 3.7, 2.0))  # socket af 7.3
    assert not a.is_inside((0, 0, 0.55)) and a.is_inside((0, 0, 0.45))  # socket floor
    assert not a.is_inside((12.5, 0, 1.0))  # screw hole
    assert a.is_inside((6.1, 0, -10.0)) and not a.is_inside((6.3, 0, -10.0))  # upper shaft
    assert a.is_inside((3.9, 0, -25.0)) and not a.is_inside((4.1, 0, -25.0))  # lower shaft
    assert a.is_inside((5.0, 0, -18.9)) and not a.is_inside((5.0, 0, -19.1))  # 45 deg step
    assert not a.is_inside((1.55, 0, -29.95)) and a.is_inside((1.65, 0, -29.95))  # magnet d 3.2
    assert not a.is_inside((0, 0, -27.95)) and a.is_inside((0, 0, -27.85))  # depth 2.1
    s = cp.hex_socket_cut(7.0, 0.15, 5.0, 5.5, 0.4)
    assert s.is_inside((0, 3.6, 7.0)) and not s.is_inside((0, 3.7, 7.0))
    assert s.is_inside((3.6 + 0.35, 0, 5.02))  # chamfered mouth
    assert s.bounding_box().max.Z == pytest.approx(10.5)


def test_bearing_seat_and_pin():
    seat = cp.bearing_seat_cut(22.0, 7.0, 0.05, 20.0, 11.0)
    assert seat.is_inside((11.0, 0, 16.0)) and not seat.is_inside((11.1, 0, 16.0))
    assert seat.is_inside((10.0, 0, 12.6)) and not seat.is_inside((10.5, 0, 12.0))  # 45 deg cone
    assert seat.bounding_box().min.Z == pytest.approx(11.0)
    pin = cp.printed_pin(3.95, 0.0, 15.0, 0.5, 5.25, 0.75)
    assert pin.is_inside((3.9, 0, 5.0)) and not pin.is_inside((4.0, 0, 5.0))
    assert not pin.is_inside((3.8, 0, 0.05))  # chamfered tip
    assert pin.is_inside((5.2, 0, 14.6)) and not pin.is_inside((5.2, 0, 14.0))  # collar
    b = cp.bearing_dummy(8, 22, 7, 0)
    assert b.volume == pytest.approx(math.pi / 4 * (22**2 - 8**2) * 7)
