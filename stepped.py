"""Parametric stepped-planet (compound) planetary reducer, ratio 10, for the GM5208-24 (build123d).

Each planet is one print with two gears on a common axis: the large gear meshes with the sun
(herringbone, module m1), the small gear meshes with the fixed ring (spur, module m2).
Ring fixed, sun input, carrier output: ratio = 1 + zr*z1 / (zs*z2) = 1 + 54*26/(12*13) = 10.

Axis = Z, z = 0 on the motor rotor's top face. Stack from -Z to +Z:
  base (screwed to the stator from below) | motor in housing_lower (tube, topped by the lip:
  the floor under the planets) | sun adapter (on the rotor) | housing_upper: ring teeth around
  the small gears and the sun shaft, large gears + sun teeth inside a thin shell with the external cap thread on
  its top band | carrier (output) with 3 printed pins | screw-on cap (lid).

Sizes are geometric lower bounds ("assembles, turns without interference, prints"), not stress
limits: size_stepped() lists them and their slack. Lewis stress, pin bending and pressures are
reference values only (reference_loads), never pass/fail.

Planet support: one 608ZZ per planet, pressed into the large gear's top end; the carrier's
printed Ø8 pin goes through its inner race and on down through a plain bore in the small gear
(a printed journal, radial clearance journal_clr). pin_kind="m8_bolt" swaps the printed pins
for M8 low-head bolts screwed into the carrier from the output face.

The housing is split at the ring teeth's bottom face (ring_lo = lip top). Both split faces are
bed faces (housing_upper prints upright, housing_lower upside down on its lip top), so the
spigot is a loose printed ring (spigot_ring) sitting half in a groove in each face, outside the
joint screw circle. 4 M3 cap screws, put in from the ring-top face through the ring band
outside the ring teeth, join the halves. The
lid is a cap screwed onto an external thread (pitch 2, 45 deg flanks, bd_warehouse) on the
housing's top band; it stops on the housing top face and a printed ratchet arm in its plate
keeps it from unscrewing. The motor (d 60) cannot pass the ring bore, so the base is separate.

Assembly order:
  1 screw sun_adapter to the rotor; 2 screw base to the stator (cable through the floor notch);
  3 lower housing_lower over motor + adapter, screw it to the base from below (4 M3x12);
    drop spigot_ring into its groove, put housing_upper on it, screw it down from inside (4 M3x12
    through the counterbores in the ring-top face; the large gears are not in yet);
  4 press a 608ZZ into the large-gear end of each planet; mesh the 3 planets radially around the
    sun (herringbone cannot slide axially), phase marks outward; lower the carrier so its 3 pins
    go through the inner races and the small-gear bores;
  5 lower this cluster: the small gears slide into the ring (turn the cluster < 1 ring pitch to
    find the mesh), the sun's hex socket slides onto the adapter hex and the sun's bottom face
    rests on the adapter's round seat (turn the rotor by hand);
  6 screw the cap on (clockwise from above) until its plate seats on the housing top; the
    arm's tooth clicks over the notches and blocks unscrewing. To open, lift the arm end.

Gear formulas (KHK Gear Technical Reference, internal gear interference) and the tooth search
are from the design study (scratchpad stepped.py); the 3-planet assembly condition
(zs*z2 + zr*z1)/3 integer is from Zhuang & Tsai, Springer, DOI 10.1007/978-3-031-45709-8_5.
"""

from __future__ import annotations

import json
import math
import sys
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
from build123d import (
    Axis,
    Box,
    Compound,
    Cone,
    Cylinder,
    Face,
    GeomType,
    Location,
    Part,
    Polyline,
    Pos,
    Rot,
    Solid,
    Wire,
    export_step,
    export_stl,
    extrude,
)

import components as cp
import planetary as pl
from planetary import half_thickness, herringbone, inv, outline, ring_half_space
from reducer import BEARINGS, MIN, SCREW_LENGTHS

M8_LENGTHS = (16, 20, 25, 30, 35, 40)  # stock M8 low-head socket screws (DIN 7984) [mm]


@dataclass
class Params:
    # gears: sun / planet large (mesh 1, module m1) and planet small / ring (mesh 2, module m2)
    zs: int = 12
    z1: int = 26
    z2: int = 13
    zr: int = 54
    n_planets: int = 3
    m1: float = 1.1
    m2: float = 1.02
    alpha_deg: float = 20.0
    helix_deg: float = 20.0  # sun-large herringbone, as planetary.py
    xs: float = 0.3  # z = 12 needs >= 0.30 against undercut
    x1: float = -0.3  # xs + x1 = 0: standard centre distance a = m1 (zs + z1) / 2 = 20.9
    x2: float = 0.4  # xr follows from coaxiality at a (about +0.39)
    # sun-large: 608 width 7 + seat cone (bearing_shoulder + fit_clr) + min_wall floor joining
    # the large gear to the small one under the cone
    b1: float = 9.25
    b2: float = 8.0  # small-ring: journal on the pin, length = pin diameter
    sun_tip_trim: float = 0.1
    small_tip_trim: float = 0.15
    flank_clr: float = 0.05
    fillet_m: float = 0.38
    ring_fillet_m: float = 0.3
    pts_per_flank: int = 8
    tip_min: float = 0.45
    # reference loads only (Lewis, Kv = 1, ring form factor of a rack 0.154)
    torque_cont: float = 0.25
    torque_peak: float = 0.4
    # motor GM5208-24 (see planetary.py for sources)
    motor_d: float = 60.0
    motor_h: float = 23.25
    shaft_od: float = 15.0
    shaft_id: float = 12.6
    rotor_n: int = 4
    rotor_pcd: float = 25.0
    rotor_phase_deg: float = 0.0
    stator_n: int = 4
    stator_pcd: float = 30.0
    stator_phase_deg: float = 0.0
    rotor_thread_depth: float = 2.7
    stator_thread_depth: float = 2.7
    motor_clr: float = 0.75
    motor_top_clr: float = 0.5  # rotor top to the housing lip
    # sun drive
    adapter_t: float = 3.3  # M3x6 into the rotor: 6 - 3.3 = 2.7 = rotor thread depth
    head_h: float = 1.3  # ESCO ultra-low socket head M3 (head d assumed 5.5 = screw_head_d)
    head_gap: float = 0.5
    hex_af: float = 7.0
    hex_engage: float = 5.0
    hex_clr: float = 0.15
    hex_chamfer: float = 0.4  # 45 deg at the sun socket mouth (first layer squish)
    sun_seat_r: float = 8.5  # adapter seat under the sun's bottom face (s0); hex only drives
    sun_shaft_r: float = 9.0
    gear_gap: float = 0.5  # small gear to lip, ring top to large gear
    sweep_clr: float = 0.5  # large tip circle to the shell bore
    # planet support: 608ZZ in the large gear, printed pin (or M8 bolt) from the carrier
    bearing: str = "608ZZ"
    bearing_shoulder: float = 1.0  # outer race rests on the seat cone this far inside D/2
    inner_land_d: float = 12.0  # 608 inner-race face OD (assumed; catalogue d1 ~ 12)
    pin_kind: str = "printed"  # or "m8_bolt"
    pin_fit_clr: float = 0.05  # pin radius d/2 - this: slip fit in the inner race
    # small-gear bore on the pin: printed hole on a printed pin, both ~0.1 off nominal; 0.05
    # (slide_clr) is for short locating fits, a 8 mm long running fit that must turn freely
    # needs about 0.2-0.3. 0.25 -> bore d 8.4 on a d 7.9 pin.
    journal_clr: float = 0.25
    pin_chamfer: float = 0.5
    seat_d: float = 10.5  # carrier collar on the inner race (inside inner_land_d)
    seat_gap: float = 0.25  # collar to the inner race; the sun sits on the adapter seat
    bolt_head_d: float = 13.0  # M8 DIN 7984
    bolt_head_h: float = 5.0
    bolt_tap_d: float = 6.8  # the bolt cuts its thread in the printed carrier
    # carrier and lid
    carrier_gap: float = 1.0
    carrier_t: float = 4.0  # >= output hole depth 8 - (out - c1) + min_wall
    carrier_wall: float = 2.0
    carrier_rim_clr: float = 1.0
    hub_step: float = 2.0
    hub_stub: float = 3.0
    n_out_holes: int = 4
    out_hole_depth: float = 8.0
    out_phase_deg: float = 45.0  # m8_bolt: use 3 holes at 60 deg, between the bolt heads
    lid_t: float = 2.7
    lid_bore_clr: float = 0.05
    lid_axial_clr: float = 0.05
    lid_rib_h: float = 1.0
    base_t: float = 3.3
    flush_cap: bool = True  # body OD = cap OD (cap flush with the body below the thread)
    # housing split at the ring bottom (lip top); loose spigot ring outside the joint screws,
    # spigot_h deep in each half, centred on the grooves' outer walls
    spigot_w: float = 1.6
    spigot_h: float = 1.5
    spigot_in_clr: float = 0.3
    n_joint: int = 4
    joint_phase_deg: float = 45.0  # between the base screws
    joint_clamp: float = 4.0
    base_clamp: float = 3.3  # = base_t: no counterbore, M3x6 still engages the stator 2.7
    wire_notch_deg: float = 45.0
    wire_notch_r0: float = 19.0
    wire_notch_w: float = 10.0
    # cap thread on the housing top band: pitch 2, 45 deg flanks (90 deg included) print without
    # support, crest/root flats 0.2. thread_clr radial on both diameters: flank normal gap
    # 0.3 sin 45 = 0.21, FDM threads need ~0.2-0.4 (slide_clr 0.05 is for short smooth bores).
    thread_pitch: float = 2.0
    thread_depth: float = 0.8
    thread_clr: float = 0.3
    thread_turns: int = 3  # cap thread length (turns + 0.5) pitch; housing ridge one turn longer
    # ratchet lock: arm in the cap plate rim, its tooth drops into notches in the housing top
    arm_deg: float = 40.0
    arm_w: float = 4.0
    arm_slot: float = 0.8
    arm_at_deg: float = 70.0  # free end of the arm (root at arm_at_deg + arm_deg)
    tooth_h: float = 0.6
    tooth_flat: float = 1.0
    notch_play: float = 0.4
    # 40 notches = 9 deg; thread_clr / pitch * 360 = 54 deg = 6 notches
    n_notches: int = 40
    n_grip: int = 48  # vertical grooves on the cap plate rim
    grip_w: float = 1.2
    grip_d: float = 0.6
    # fasteners and FDM
    n_base_screws: int = 4
    m3_clear_d: float = 3.4
    m3_tap_d: float = 2.5
    screw_d: float = 3.0
    screw_engage: float = 6.0
    hex_hole_extra: float = 1.0
    screw_head_d: float = 5.5
    joint_head_h: float = 3.0  # ISO 4762 M3 socket head
    fit_clr: float = 0.05
    slide_clr: float = 0.05
    min_wall: float = 1.2

    # --- derived: gears ---
    @property
    def alpha(self) -> float:
        return math.radians(self.alpha_deg)

    @property
    def ratio(self) -> float:
        return 1 + self.zr * self.z1 / (self.zs * self.z2)

    @property
    def a(self) -> float:
        return self.m1 * (self.zs + self.z1) / 2

    @property
    def aw2(self) -> float:
        """Working pressure angle of the ring mesh at centre distance a."""
        return math.acos(self.m2 * (self.zr - self.z2) / 2 * math.cos(self.alpha) / self.a)

    @property
    def xr(self) -> float:
        """Ring shift for coaxiality: inv aw = inv a + 2 tan a (xr - x2)/(zr - z2)."""
        d = (inv(self.aw2) - inv(self.alpha)) * (self.zr - self.z2) / (2 * math.tan(self.alpha))
        return self.x2 + d

    def mesh(self, k: int) -> pl.Params:
        """planetary.Params carrying mesh k's module/shifts, for the planetary gear helpers."""
        if k == 1:
            return pl.Params(zs=self.zs, zp=self.z1, m=self.m1, xs=self.xs, xp=self.x1,
                             face_w=self.b1, flank_clr=self.flank_clr, fillet_m=self.fillet_m,
                             pts_per_flank=self.pts_per_flank, helix_deg=self.helix_deg)  # fmt: skip
        return pl.Params(zp=self.z2, zr=self.zr, m=self.m2, xp=self.xr, flank_clr=self.flank_clr,
                         fillet_m=self.fillet_m, pts_per_flank=self.pts_per_flank)  # fmt: skip

    @property
    def r_tip_sun(self) -> float:
        return self.m1 * (self.zs / 2 + 1 + self.xs) - self.sun_tip_trim

    @property
    def r_tip1(self) -> float:
        return self.m1 * (self.z1 / 2 + 1 + self.x1)

    @property
    def r_tip2(self) -> float:
        return self.m2 * (self.z2 / 2 + 1 + self.x2) - self.small_tip_trim

    def r_root(self, m: float, z: int, x: float) -> float:
        return m * (z / 2 - 1.25 + x)

    @property
    def r_ring_tip(self) -> float:
        return self.m2 * (self.zr / 2 - 1 + self.xr)

    @property
    def r_ring_root(self) -> float:
        return self.m2 * (self.zr / 2 + 1.25 + self.xr)

    @property
    def R_sweep(self) -> float:  # large-gear swing circle + sweep_clr
        return self.a + self.r_tip1 + self.sweep_clr

    # --- derived: planet support ---
    @property
    def pin_r(self) -> float:
        return BEARINGS[self.bearing][0] / 2 - self.pin_fit_clr

    @property
    def journal_r(self) -> float:  # small-gear bore
        return self.pin_r + self.journal_clr

    @property
    def seat_r(self) -> float:  # 608 seat in the large gear
        return BEARINGS[self.bearing][1] / 2 + self.fit_clr

    # --- derived: housing radii (each the smallest value its constraints allow) ---
    @property
    def R_motor_bore(self) -> float:
        return self.motor_d / 2 + self.motor_clr

    @property
    def hex_circ(self) -> float:  # circumradius of an M3 hex hole
        return self.screw_d / 2 / math.cos(math.pi / 6)

    @property
    def cbore_r(self) -> float:
        return self.screw_head_d / 2 + self.slide_clr

    @property
    def joint_r(self) -> float:
        """Joint screws: counterbore min_wall outside the ring root, hex hole outside the motor."""
        return max(self.r_ring_root + self.min_wall + self.cbore_r,
                   self.R_motor_bore + self.min_wall + self.hex_circ)  # fmt: skip

    @property
    def R_shell(self) -> float:
        """Shell bore above the ring: the large-gear swing and the joint counterbore columns."""
        return max(self.R_sweep, self.joint_r + self.cbore_r)

    @property
    def R_thread_root(self) -> float:  # min_wall behind the external thread
        return self.R_shell + self.min_wall

    @property
    def R_crest(self) -> float:
        return self.R_thread_root + self.thread_depth

    @property
    def R_body(self) -> float:
        r = max(self.R_crest, self.joint_r + self.cbore_r + self.min_wall,
                self.joint_r + self.hex_circ + self.min_wall, self.R_motor_bore + 2.0,
                self.spigot_r[1] + self.slide_clr + self.min_wall)  # fmt: skip
        return max(r, self.R_cap) if self.flush_cap else r

    @property
    def R_cap(self) -> float:
        return cp.cap_outer_r(self.R_thread_root, self.thread_depth, self.thread_clr,
                              self.min_wall)  # fmt: skip

    @property
    def base_screw_r(self) -> float:
        """Midway between the motor-bore bound and the outside bound of the base screw holes."""
        lo = self.R_motor_bore + self.min_wall + self.hex_circ
        hi = self.R_body - self.min_wall - max(self.hex_circ, self.cbore_r)
        return (lo + hi) / 2

    @property
    def spigot_r(self) -> tuple[float, float]:
        """(inner, outer) radius of the spigot ring: groove min_wall outside the joint holes."""
        r0 = (self.joint_r + max(self.m3_clear_d / 2, self.hex_circ) + self.min_wall
              + self.spigot_in_clr)  # fmt: skip
        return r0, r0 + self.spigot_w

    @property
    def adapter_r(self) -> float:
        return self.rotor_pcd / 2 + self.screw_head_d / 2 + 1.1 * self.min_wall

    @property
    def lip_r(self) -> float:  # housing lip bore: the rotating adapter passes through it
        return self.adapter_r + 1.0

    @property
    def hex_sock(self) -> float:
        return (self.hex_af / 2 + self.hex_clr) / math.cos(math.pi / 6)

    @property
    def bolt_cbore_r(self) -> float:
        return self.bolt_head_d / 2 + self.slide_clr

    @property
    def carrier_r(self) -> float:
        r = self.seat_d / 2 + self.carrier_wall
        if self.pin_kind == "m8_bolt":  # bolt-head counterbore inside the hub
            r = max(r, self.bolt_cbore_r + self.min_wall + self.hub_step)
        return self.a + r

    @property
    def hub_r(self) -> float:
        return self.carrier_r - self.hub_step

    @property
    def cavity_r(self) -> float:
        return self.carrier_r + self.carrier_rim_clr

    @property
    def zs_(self) -> dict[str, float]:
        """Key heights: s0/s1 small gear (ring teeth ring_lo..s1), L0/L1 large gear."""
        s0 = max(
            self.adapter_t + self.head_h + self.head_gap,  # sun on the adapter
            self.motor_top_clr + self.min_wall + self.gear_gap,
        )  # planets over the lip
        s1 = s0 + self.b2
        L0 = s1 + self.gear_gap
        L1 = L0 + self.b1
        c0 = L1 + self.carrier_gap
        c1 = c0 + self.carrier_t
        lid = c1 + self.lid_axial_clr + self.lid_rib_h
        base_top = -self.motor_h
        return {"s0": s0, "ring_lo": s0 - self.gear_gap, "s1": s1, "L0": L0, "L1": L1,
                "Lmid": (L0 + L1) / 2, "c0": c0, "c1": c1, "lid": lid, "top": lid + self.lid_t,
                "out": lid + self.lid_t + self.hub_stub,
                "thread0": lid - (self.thread_turns + 1) * self.thread_pitch,
                "base_top": base_top, "bottom": base_top - self.base_t,
                "split": s0 - self.gear_gap,
                "joint_seat": s0 - self.gear_gap + self.joint_clamp}  # fmt: skip

    @property
    def height(self) -> float:
        return self.zs_["out"] - self.zs_["bottom"]

    @property
    def lip_top(self) -> float:  # housing lip under the small gears
        return self.zs_["s0"] - self.gear_gap

    @property
    def bolt_len(self) -> int:
        """Longest stock M8 that ends at or above s0 (m8_bolt variant)."""
        need = self.zs_["out"] - self.bolt_head_h - self.zs_["s0"]
        return max(L for L in M8_LENGTHS if L <= need)

    def screw_len(self, clamp: float) -> int:
        return next(L for L in SCREW_LENGTHS if L >= clamp + self.screw_engage)

    def motor_screw_len(self, clamp: float, depth: float) -> int:
        return max((L for L in SCREW_LENGTHS if L <= clamp + depth), default=0)


# ---------------------------------------------------------------- checks
def xmin(p: Params, z: int) -> float:
    return 1 - z * math.sin(p.alpha) ** 2 / 2


def tip_land(p: Params, m: float, z: int, x: float, ra: float) -> float:
    """External tooth chord at radius ra, flank_clr removed (planetary.tip_land, any module)."""
    q = pl.Params(m=m, flank_clr=p.flank_clr)
    return pl.tip_land(q, z, x, ra)


def ring_tip_land(p: Params) -> float:
    m, z, ra = p.m2, p.zr, p.r_ring_tip
    rb = m * z / 2 * math.cos(p.alpha)
    e = ring_half_space(p.mesh(2))
    return 2 * ra * (math.pi / z - (e + inv(p.alpha) - inv(math.acos(rb / ra))))


def internal_interference(p: Params) -> dict[str, bool]:
    """KHK Gear Technical Reference, internal gears: True = interferes."""
    m, z1, z2, a, aw = p.m2, p.z2, p.zr, p.a, p.aw2
    ra1, ra2 = p.r_tip2, p.r_ring_tip
    rb1, rb2 = m * z1 / 2 * math.cos(p.alpha), m * z2 / 2 * math.cos(p.alpha)
    aa1, aa2 = math.acos(rb1 / ra1), math.acos(rb2 / ra2)
    involute = z1 / z2 < 1 - math.tan(aa2) / math.tan(aw)
    c1 = (ra2**2 - ra1**2 - a**2) / (2 * a * ra1)
    c2 = (a**2 + ra2**2 - ra1**2) / (2 * a * ra2)
    trochoid = False
    if abs(c1) <= 1 and abs(c2) <= 1:
        th1 = math.acos(c1) + inv(aa1) - inv(aw)
        trochoid = th1 * z1 / z2 + inv(aw) - inv(aa2) - math.acos(c2) < 0
    return {"involute": involute, "trochoid": trochoid}


def planet_tip_gap(p: Params) -> float:
    return 2 * p.a * math.sin(math.pi / p.n_planets) - 2 * p.r_tip1


def lewis_stress(p: Params, torque_in: float) -> dict[str, float]:
    """Lewis root stress [MPa], Kv = 1 (reference only). Ring side Ft2 = Ft1 r1/r2."""
    ft1 = torque_in / (p.m1 * p.zs / 2e3) / p.n_planets
    ft2 = ft1 * (p.m1 * p.z1) / (p.m2 * p.z2)

    def y(z: int) -> float:
        return 0.154 - 0.912 / z

    return {"sun": ft1 / (p.b1 * y(p.zs) * math.pi * p.m1),
            "planet_large": ft1 / (p.b1 * y(p.z1) * math.pi * p.m1),
            "planet_small": ft2 / (p.b2 * y(p.z2) * math.pi * p.m2),
            "ring": ft2 / (p.b2 * 0.154 * math.pi * p.m2)}  # fmt: skip


def reference_loads(p: Params, torque_in: float) -> dict[str, float]:
    """Rough planet-support loads at input torque_in [N m] (reference values, not checks).

    Per planet the sun force Ft1 acts at the large gear's mid-plane and the ring force Ft2 at
    the small gear's; both push the planet the same way (their sum drives the carrier). The
    planet rests on the 608 (centre zb) and the journal (centre zj = small mid-plane): moment
    balance about zj gives the 608 load Rb, the rest Rj goes through the journal. The pin is a
    cantilever from the carrier underside c0 (collar ignored): M = Rb (c0 - zb) + Rj (c0 - zj).
    Pressures are projected (load / (d x length)).
    """
    z = p.zs_
    B = BEARINGS[p.bearing][2]
    d = 2 * p.pin_r
    ft1 = torque_in / (p.m1 * p.zs / 2e3) / p.n_planets
    ft2 = ft1 * (p.m1 * p.z1) / (p.m2 * p.z2)
    zb, zj = z["L1"] - B / 2, (z["s0"] + z["s1"]) / 2
    rb = ft1 * (z["Lmid"] - zj) / (zb - zj)
    rj = ft1 + ft2 - rb
    m = rb * (z["c0"] - zb) + rj * (z["c0"] - zj)
    return {"Ft1 [N]": ft1, "Ft2 [N]": ft2, "608 load [N]": rb, "journal load [N]": rj,
            "pin root moment [N mm]": m, "pin bending [MPa]": 32 * m / (math.pi * d**3),
            "journal pressure [MPa]": rj / (d * p.b2),
            "carrier pin-root pressure [MPa]": (ft1 + ft2) / (d * p.carrier_t)}  # fmt: skip


def design_warnings(p: Params) -> list[str]:
    """Geometric problems only (stress is reported by reference_loads, never checked)."""
    w = []
    z = p.zs_
    B = BEARINGS[p.bearing][2]
    if (p.zs * p.z2 + p.zr * p.z1) % p.n_planets:
        w.append("assembly: (zs z2 + zr z1)/N not integer")
    if abs(p.xs + p.x1) > 1e-9:
        w.append("xs + x1 != 0: model assumes the standard sun-large centre distance")
    if not 0.80 < math.cos(p.aw2) < 0.99:
        w.append(f"ring mesh working pressure angle {math.degrees(p.aw2):.1f} deg out of range")
    if planet_tip_gap(p) < p.min_wall:
        w.append(f"planet tips {planet_tip_gap(p):.2f} mm apart")
    for name, m, zz, x, ra in (("sun", p.m1, p.zs, p.xs, p.r_tip_sun),
                               ("planet_large", p.m1, p.z1, p.x1, p.r_tip1),
                               ("planet_small", p.m2, p.z2, p.x2, p.r_tip2)):  # fmt: skip
        if x < xmin(p, zz) - 1e-9:
            w.append(f"undercut: {name} z={zz} x={x}")
        if tip_land(p, m, zz, x, ra) < p.tip_min:
            w.append(f"{name} tip land {tip_land(p, m, zz, x, ra):.2f} < {p.tip_min}")
    if ring_tip_land(p) < p.tip_min:
        w.append(f"ring tip land {ring_tip_land(p):.2f} < {p.tip_min}")
    if p.r_ring_tip <= p.m2 * p.zr / 2 * math.cos(p.alpha):
        w.append("ring tip inside its base circle")
    w += [f"internal {k} interference" for k, v in internal_interference(p).items() if v]
    if p.sun_seat_r > min(p.sun_shaft_r, p.rotor_pcd / 2 - p.screw_head_d / 2 - p.slide_clr):
        w.append("adapter sun seat wider than the sun shaft or reaching the rotor screw heads")
    if p.a - p.r_tip2 - p.sun_shaft_r < 1.0:
        w.append("sun shaft within 1 mm of the small gear tips")
    if p.b1 - B < p.bearing_shoulder + p.fit_clr + p.min_wall - 1e-9:
        w.append("large gear too thin for the 608 and its seat cone")
    if p.seat_d > p.inner_land_d or p.seat_d / 2 <= p.pin_r:
        w.append("carrier collar not on the inner-race land")
    if p.lip_top < p.motor_top_clr + p.min_wall - 1e-9:
        w.append("housing lip thinner than min_wall")
    if z["split"] < z["ring_lo"] - 1e-9:
        w.append("housing split below the ring teeth: the lip would be in housing_upper")
    if p.spigot_r[1] + p.slide_clr + p.min_wall > p.R_body + 1e-9:
        w.append("spigot groove within min_wall of the outside")
    if z["joint_seat"] + p.joint_head_h > z["s1"] + 1e-9:
        w.append("joint screw heads stand above the ring top")
    if p.R_body > p.R_cap + 1e-9:
        w.append("housing body wider than the cap")
    if p.screw_len(p.joint_clamp) - p.joint_clamp < 4:
        w.append("joint screws engage < 4 mm in the tube")
    if p.R_motor_bore < p.r_ring_tip:
        w.append("motor would pass the ring bore: base could be part of the housing")
    if z["thread0"] < z["s1"]:
        w.append("cap thread zone reaches below the ring top")
    for face, clamp, depth in (("rotor", p.adapter_t, p.rotor_thread_depth),
                               ("stator", p.base_clamp, p.stator_thread_depth)):  # fmt: skip
        if p.motor_screw_len(clamp, depth) - clamp < 2:
            w.append(f"{face} screws engage < 2 mm (no stock length fits)")
    if p.R_cap - p.arm_w > p.R_shell + 0.1:
        w.append("ratchet arm does not reach over the housing top face")
    if p.pin_kind == "m8_bolt":
        if z["out"] - p.bolt_head_h - p.bolt_len > z["s1"] - 1.0:
            w.append(f"M8x{p.bolt_len} ends above the small gear's journal")
        ro = p.hub_r / 2
        for k in range(p.n_out_holes):
            q = ro * np.exp(1j * math.radians(p.out_phase_deg + 360 * k / p.n_out_holes))
            for j in range(p.n_planets):
                c = p.a * np.exp(2j * math.pi * j / p.n_planets)
                if abs(q - c) < p.bolt_cbore_r + p.m3_tap_d / 2 + 0.5:
                    w.append("output hole runs into an M8 head counterbore")
    elif p.pin_kind != "printed":
        w.append(f"unknown pin_kind {p.pin_kind}")
    return w


def wall_thicknesses(p: Params) -> dict[str, float]:
    z = p.zs_
    return {
        "sun shaft around hex socket": p.sun_shaft_r - p.hex_sock,
        "planet small root to journal bore": p.r_root(p.m2, p.z2, p.x2) - p.journal_r,
        "planet large root to 608 seat": p.r_root(p.m1, p.z1, p.x1) - p.seat_r,
        "ring root to outside": p.R_body - p.r_ring_root,
        "shell behind the cap thread": p.R_thread_root - p.R_shell,
        "housing around motor": p.R_body - p.R_motor_bore,
        "lower: lip (axial)": p.lip_top - p.motor_top_clr,
        "upper: joint counterbore to ring root": p.joint_r - p.cbore_r - p.r_ring_root,
        "upper: joint counterbore to outside": p.R_body - p.joint_r - p.cbore_r,
        "spigot groove to joint hole": p.spigot_r[0] - p.spigot_in_clr - p.joint_r
        - max(p.m3_clear_d / 2, p.hex_circ),
        "spigot groove to outside": p.R_body - p.spigot_r[1] - p.slide_clr,
        "upper: under the joint heads": p.joint_clamp,
        "lower: joint hole to motor bore": p.joint_r - p.hex_circ - p.R_motor_bore,
        "lower: joint hole to outside": p.R_body - p.joint_r - p.hex_circ,
        "housing: base screw hole to motor bore": p.base_screw_r - p.hex_circ - p.R_motor_bore,
        "housing: base screw hole to outside": p.R_body - p.base_screw_r - p.hex_circ,
        "base: screw head to outside": p.R_body - p.base_screw_r - p.cbore_r,
        "carrier plate around collars": p.carrier_r - p.a - p.seat_d / 2,
        "cap plate": p.lid_t,
        "cap skirt": p.min_wall,
        "adapter around screw holes": p.adapter_r - p.rotor_pcd / 2 - p.m3_clear_d / 2,
        "adapter boss around rotor screw heads": p.rotor_pcd / 2 - p.screw_head_d / 2 - p.hex_sock,
        "base under counterbores": p.base_clamp,
        "base: centre hole to stator screw": p.stator_pcd / 2 - p.m3_clear_d / 2
        - (p.shaft_od / 2 + 0.5),
        "base: wire notch to stator counterbore": p.stator_pcd / 2 * math.sin(math.radians(
            min(p.wire_notch_deg % (360 / p.stator_n), -p.wire_notch_deg % (360 / p.stator_n))))
        - p.wire_notch_w / 2 - p.cbore_r,
        "hub around output holes": p.hub_r / 2 - p.m3_tap_d / 2,
        "carrier under the output holes": z["out"] - p.out_hole_depth - z["c0"],
        "planet tips": planet_tip_gap(p),
    }  # fmt: skip


def thin_walls(p: Params) -> dict[str, float]:
    return {k: v for k, v in wall_thicknesses(p).items() if v < p.min_wall - 1e-9}


def size_stepped(p: Params) -> dict:
    """Outer size and the geometric lower bounds that set it (slack >= 0 means satisfied).

    Replaces the sizing study's size_stepped(zs, z1, z2, zr, m1, m2, xs, x2, bearing_in="large"):
    the same constraints, now read from Params so the summary and the model agree.
    """
    z = p.zs_
    B = BEARINGS[p.bearing][2]
    cand = {"cap thread crest (shell + min_wall + depth)": p.R_crest,
            "joint counterbore + min_wall": p.joint_r + p.cbore_r + p.min_wall,
            "joint hex hole + min_wall": p.joint_r + p.hex_circ + p.min_wall,
            "motor bore + 2": p.R_motor_bore + 2.0,
            "spigot groove + min_wall": p.spigot_r[1] + p.slide_clr + p.min_wall}  # fmt: skip
    if p.flush_cap:
        cand["cap OD (flush_cap)"] = p.R_cap
    shell = {"large-gear swing + sweep_clr": p.R_sweep,
             "joint counterbore column": p.joint_r + p.cbore_r}  # fmt: skip
    joint = {"ring root + min_wall + counterbore": p.r_ring_root + p.min_wall + p.cbore_r,
             "motor bore + min_wall + hex hole": p.R_motor_bore + p.min_wall + p.hex_circ}  # fmt: skip
    s0 = {"adapter + rotor screw head + gap": p.adapter_t + p.head_h + p.head_gap,
          "motor top clr + lip min_wall + gear_gap": p.motor_top_clr + p.min_wall + p.gear_gap}  # fmt: skip
    height = {"base": p.base_t, "motor": p.motor_h, "s0 (sun floor)": z["s0"],
              "b2 small gear (journal = pin d)": p.b2, "gear_gap": p.gear_gap,
              "b1 large gear (608 B + seat cone + floor)": p.b1, "carrier_gap": p.carrier_gap,
              "carrier_t": p.carrier_t, "lid rib + axial clr": p.lid_rib_h + p.lid_axial_clr,
              "lid_t": p.lid_t, "hub_stub": p.hub_stub}  # fmt: skip
    slack = {
        "sun shaft - hex socket - min_wall": p.sun_shaft_r - p.hex_sock - p.min_wall,
        "small root - journal bore - min_wall": p.r_root(p.m2, p.z2, p.x2) - p.journal_r
        - p.min_wall,
        "large root - 608 seat - min_wall": p.r_root(p.m1, p.z1, p.x1) - p.seat_r - p.min_wall,
        "planet tips - min_wall": planet_tip_gap(p) - p.min_wall,
        "b1 - (B + seat cone + min_wall floor)": p.b1 - B - p.bearing_shoulder - p.fit_clr
        - p.min_wall,
        "lower joint hole - motor bore - min_wall": p.joint_r - p.hex_circ - p.R_motor_bore
        - p.min_wall,
        "upper joint counterbore - ring root - min_wall": p.joint_r - p.cbore_r - p.r_ring_root
        - p.min_wall,
        "carrier under output holes - min_wall": z["out"] - p.out_hole_depth - z["c0"]
        - p.min_wall,
    }  # fmt: skip
    return {"D_body": 2 * p.R_body, "D_max": 2 * p.R_cap, "height": p.height,
            "R_body from": cand, "R_shell from": shell, "joint_r from": joint, "s0 from": s0,
            "height terms": height, "slack": slack}  # fmt: skip


# ---------------------------------------------------------------- kinematics
def _sun_phase(p: Params, th: float) -> float:
    """Planet turn at carrier angle th whose large teeth fill the sun spaces (planetary.py)."""
    u = (p.zs * th / (2 * math.pi)) % 1
    return th + math.pi - 2 * math.pi * (0.5 - u) / p.z1


def ring_phase(p: Params) -> float:
    """Turn of the ring cutter so that the ring teeth fill planet 0's small-gear spaces."""
    v = p.z2 * (0 - _sun_phase(p, 0.0)) / (2 * math.pi)
    return -2 * math.pi * (v - 0.5) / p.zr + math.pi / p.zr


def planet_phase(p: Params, k: int) -> tuple[float, float]:
    """(carrier angle, own turn) of planet k at input 0, meshing with both sun and ring.

    Tooth 0 of the large and of the small gear lie on one radial line (phase mark), so all planets
    are the same print. The large gear may sit any whole large pitch further on; the assembly
    condition guarantees one of these also fits the ring.
    """
    th = 2 * math.pi * k / p.n_planets
    # small-gear turn the ring needs at th (ring_phase solved for the planet turn)
    need = th - 2 * math.pi * (1 - (ring_phase(p) - th) * p.zr / (2 * math.pi)) / p.z2
    pitch2 = 2 * math.pi / p.z2

    def miss(j: int) -> float:
        e = (_sun_phase(p, th) + 2 * math.pi * j / p.z1 - need) % pitch2
        return min(e, pitch2 - e)

    j = min(range(p.z1), key=miss)
    if miss(j) > 1e-9:
        raise ValueError(f"planet {k} cannot mesh with sun and ring: assembly condition")
    return th, _sun_phase(p, th) + 2 * math.pi * j / p.z1


def planet_states(p: Params, psi: float) -> list[tuple[complex, float]]:
    pc = psi / p.ratio
    out = []
    for k in range(p.n_planets):
        th, phi = planet_phase(p, k)
        out.append((p.a * np.exp(1j * (th + pc)), phi + pc - (psi - pc) * p.zs / p.z1))
    return out


# ---------------------------------------------------------------- solids
def _spur(pts, z0: float, z1: float) -> Solid:
    return Pos(0, 0, z0) * extrude(Face(Wire(Polyline(*pts, close=True).edges())), z1 - z0)


def _herring(p: Params, pts, r_pitch: float, hand: int) -> Solid:
    """planetary.herringbone over L0..L1 (it puts the apex at its own mid-plane: shift)."""
    q = p.mesh(1)
    dz = p.zs_["Lmid"] - q.zs_["mid"]
    return Pos(0, 0, dz) * herringbone(q, pts, r_pitch, hand, p.zs_["L0"] - dz, p.zs_["L1"] - dz)


def tooth_r(p: Params) -> tuple[float, float]:
    """Radial extent of the ratchet tooth: over the housing top face, inside the thread root."""
    return p.R_shell + 0.1, p.R_thread_root - 0.1


def tooth_deg(p: Params) -> float:
    """Angle of the tooth's vertical (loosening-side) face, cap at its assembled turn."""
    r = sum(tooth_r(p)) / 2
    return p.arm_at_deg + math.degrees((p.tooth_flat + p.tooth_h + 0.3) / r)


def lid_part(p: Params, turn_deg: float = 0.0, lift: float = 0.0) -> Part:
    """Screw-on cap. Seated, its plate underside rests on the housing top (z["lid"]), which fixes
    the rib-to-carrier gap lid_axial_clr whatever the thread does. turn_deg/lift move it."""
    z = p.zs_
    lid = cp.screw_cap(p.R_thread_root, p.thread_depth, p.thread_pitch, p.thread_clr,
                       p.min_wall, z["lid"], p.lid_t, p.thread_turns, p.cavity_r,
                       p.lid_rib_h)  # fmt: skip
    g = p.lid_bore_clr
    bore = p.hub_r + g
    lid -= Pos(0, 0, z["c1"]) * Cone(p.carrier_r + g * math.sqrt(2), bore,
                                     p.hub_step + g * (math.sqrt(2) - 1), align=MIN)  # fmt: skip
    lid -= Pos(0, 0, z["c1"]) * Cylinder(bore, z["top"] - z["c1"], align=MIN)
    skirt_h = (p.thread_turns + 0.5) * p.thread_pitch
    lid = cp.ratchet_arm(lid, p.R_cap, p.arm_w, p.arm_slot, p.arm_at_deg, p.arm_deg, z["lid"],
                         p.lid_t, skirt_h, tooth_deg(p), tooth_r(p), p.tooth_h,
                         p.tooth_flat)  # fmt: skip
    lid -= cp.grip_grooves_cut(p.R_cap, z["lid"], p.lid_t, p.n_grip, p.grip_w, p.grip_d)
    if turn_deg or lift:
        lid = Pos(0, 0, lift) * lid.rotate(Axis.Z, turn_deg)
    return lid


def make_parts(p: Params) -> dict[str, Part]:
    """Parts at input 0, assembled. Names starting with '_' are bought."""
    z = p.zs_
    d, D, B = BEARINGS[p.bearing]
    q1, q2 = p.mesh(1), p.mesh(2)
    parts: dict[str, Part] = {}

    motor = Pos(0, 0, -p.motor_h) * Cylinder(p.motor_d / 2, p.motor_h, align=MIN)
    motor += Pos(0, 0, -p.motor_h) * Cylinder(p.shaft_od / 2, p.motor_h, align=MIN)
    parts["_motor_GM5208-24"] = motor - Pos(0, 0, -p.motor_h) * Cylinder(
        p.shaft_id / 2, p.motor_h, align=MIN
    )

    parts["sun_adapter"] = cp.sun_adapter(p.adapter_r, p.adapter_t, p.sun_seat_r, z["s0"],
                                          p.hex_af, z["s0"] + p.hex_engage, p.rotor_n,
                                          p.rotor_pcd, p.m3_clear_d, p.rotor_phase_deg)  # fmt: skip

    # sun: herringbone teeth at the large-gear level on a shaft down to the adapter hex
    sun = _herring(p, outline(q1, p.zs, half_thickness(q1, p.zs, p.xs), p.r_tip_sun,
                              p.r_root(p.m1, p.zs, p.xs)), q1.r(p.zs), 1)  # fmt: skip
    zt = z["L0"] - p.gear_gap
    # the L0 tooth section continued straight down to the shaft: no overhanging tooth bottoms
    sun += extrude(sun.faces().sort_by(Axis.Z)[0], p.gear_gap, (0, 0, -1))
    sun += Pos(0, 0, z["s0"]) * Cylinder(p.sun_shaft_r, zt - z["s0"], align=MIN)
    parts["sun"] = sun - cp.hex_socket_cut(p.hex_af, p.hex_clr, z["s0"],
                                           p.hex_engage + p.gear_gap, p.hex_chamfer)  # fmt: skip

    # planet: spur small gear s0..L0 + herringbone large gear L0..L1, tooth 0 of both on +X;
    # 608 in the large end (seat cone, min_wall floor over L0), journal bore through both
    small = _spur(outline(q2, p.z2, half_thickness(q2, p.z2, p.x2), p.r_tip2,
                          p.r_root(p.m2, p.z2, p.x2)), z["s0"], z["L0"])  # fmt: skip
    large = _herring(p, outline(q1, p.z1, half_thickness(q1, p.z1, p.x1), p.r_tip1,
                                p.r_root(p.m1, p.z1, p.x1)), q1.r(p.z1), -1)  # fmt: skip
    planet = small.fuse(large).clean()
    planet -= cp.bearing_seat_cut(D, B, p.fit_clr, z["L1"], z["L0"] + p.min_wall)
    planet -= Pos(0, 0, z["s0"]) * Cylinder(p.journal_r, z["L1"] - z["s0"], align=MIN)
    mark_r = (p.seat_r + p.r_root(p.m1, p.z1, p.x1)) / 2  # phase mark on the large top face
    planet -= Pos(mark_r, 0, z["L1"] - 0.4) * Cylinder(0.5, 0.4, align=MIN)

    carrier = Pos(0, 0, z["c0"]) * Cylinder(p.carrier_r, p.carrier_t, align=MIN)
    carrier += Pos(0, 0, z["c1"]) * Cone(p.carrier_r, p.hub_r, p.hub_step, align=MIN)
    zh = z["c1"] + p.hub_step
    carrier += Pos(0, 0, zh) * Cylinder(p.hub_r, z["out"] - zh, align=MIN)
    carrier -= cp.ring_holes(p.n_out_holes, p.hub_r / 2, p.m3_tap_d, z["out"] - p.out_hole_depth,
                             p.out_hole_depth, p.out_phase_deg)  # fmt: skip
    zc = z["L1"] + p.seat_gap  # collar bottom, seat_gap over the inner race
    if p.pin_kind == "printed":
        pin = cp.printed_pin(p.pin_r, z["s0"], z["c0"], p.pin_chamfer, p.seat_d / 2,
                             z["c0"] - zc)  # fmt: skip
    else:
        pin = Pos(0, 0, zc) * Cylinder(p.seat_d / 2, z["c0"] - zc, align=MIN)
        zb = z["out"] - p.bolt_head_h
        bolt = Pos(0, 0, zb - p.bolt_len) * Cylinder(d / 2, p.bolt_len, align=MIN)
        bolt += Pos(0, 0, zb) * Cylinder(p.bolt_head_d / 2, p.bolt_head_h, align=MIN)
    for k, (c, turn) in enumerate(planet_states(p, 0.0)):
        x, y = c.real, c.imag
        parts[f"planet_{k}"] = Pos(x, y, 0) * planet.rotate(Axis.Z, math.degrees(turn))
        parts[f"_bearing_{p.bearing}_{k}"] = Pos(x, y, 0) * cp.bearing_dummy(d, D, B, z["L1"] - B)
        carrier += Pos(x, y, 0) * pin
        if p.pin_kind == "m8_bolt":
            carrier -= Pos(x, y, zc) * Cylinder(p.bolt_tap_d / 2, z["out"] - zc, align=MIN)
            carrier -= Pos(x, y, zb) * Cylinder(p.bolt_cbore_r, p.bolt_head_h, align=MIN)
            parts[f"_bolt_M8x{p.bolt_len}_{k}"] = Pos(x, y, 0) * bolt
    parts["carrier_output"] = carrier

    # base: under the stator, screwed to it and to the housing from below
    base = Pos(0, 0, z["bottom"]) * Cylinder(p.R_body, p.base_t, align=MIN)
    base -= Pos(0, 0, z["bottom"]) * Cylinder(p.shaft_od / 2 + 0.5, p.base_t, align=MIN)
    cb = p.base_t - p.base_clamp
    for n, r in ((p.stator_n, p.stator_pcd / 2), (p.n_base_screws, p.base_screw_r)):
        base -= cp.ring_holes(n, r, p.m3_clear_d, z["bottom"], p.base_t, p.stator_phase_deg)
        if cb > 0:  # heads stand proud of the base when base_t == base_clamp
            base -= cp.ring_holes(n, r, 2 * p.cbore_r, z["bottom"], cb, p.stator_phase_deg)
    base -= (Rot(0, 0, p.stator_phase_deg + p.wire_notch_deg)
             * Pos(p.wire_notch_r0, -p.wire_notch_w / 2, z["bottom"])
             * Box(p.R_body - p.wire_notch_r0 + 1, p.wire_notch_w, p.base_t, align=MIN))  # fmt: skip
    parts["base"] = base

    # housing: motor tube, lip (in the lower), spur ring, thin shell, external cap thread on the top band
    zt0 = z["thread0"]
    h = Pos(0, 0, z["base_top"]) * Cylinder(p.R_body, zt0 - z["base_top"], align=MIN)
    h += Pos(0, 0, zt0) * Cylinder(p.R_thread_root + cp.ROOT_BURY, z["lid"] - zt0, align=MIN)
    h += cp.housing_thread(p.R_thread_root, p.thread_depth, p.thread_pitch, z["lid"],
                           p.thread_turns + 1)  # fmt: skip
    h -= Pos(0, 0, z["base_top"]) * Cylinder(p.R_motor_bore, p.motor_h + p.motor_top_clr,
                                             align=MIN)  # fmt: skip
    h -= Cylinder(p.lip_r, z["s1"], align=MIN)
    qr = replace(q2, fillet_m=p.ring_fillet_m)
    ring = _spur(outline(qr, p.zr, ring_half_space(qr), p.r_ring_root, p.r_ring_tip, cutter=True),
                 z["ring_lo"], z["s1"])  # fmt: skip
    h -= ring.rotate(Axis.Z, math.degrees(ring_phase(p)))
    h -= Pos(0, 0, z["s1"]) * Cylinder(p.R_shell, z["lid"] - z["s1"], align=MIN)
    h -= cp.ratchet_notches_cut(p.R_shell - 1, p.R_crest + 1, z["lid"], p.n_notches, tooth_deg(p),
                                sum(tooth_r(p)) / 2, p.tooth_h, p.tooth_flat,
                                p.notch_play)  # fmt: skip
    hd = p.screw_len(p.base_clamp) - p.base_clamp + p.hex_hole_extra
    for k in range(p.n_base_screws):
        ph = math.radians(p.stator_phase_deg) + 2 * math.pi * k / p.n_base_screws
        h -= cp.hex_prism(p.screw_d, z["base_top"], hd, p.base_screw_r * math.cos(ph),
                          p.base_screw_r * math.sin(ph))  # fmt: skip
    upper, lower, ring = cp.split_housing(
        h, z["split"], p.R_body, p.spigot_r[0], p.spigot_w, p.spigot_h, p.spigot_in_clr,
        p.slide_clr, p.n_joint, p.joint_r, p.joint_phase_deg, p.joint_clamp, z["s1"], p.screw_d,
        p.m3_clear_d, p.screw_head_d, p.screw_len(p.joint_clamp), p.hex_hole_extra)  # fmt: skip
    parts["housing_upper"] = upper
    parts["housing_lower"] = lower
    parts["spigot_ring"] = ring
    parts["lid"] = lid_part(p)
    return parts


def pose(p: Params, home: dict[str, Part], psi: float) -> dict[str, Part]:
    """Sun at angle psi [rad] (ring fixed, carrier psi/ratio)."""
    rot = Location((0, 0, 0), (0, 0, math.degrees(psi)))
    rc = Location((0, 0, 0), (0, 0, math.degrees(psi / p.ratio)))
    move = {"sun": rot, "sun_adapter": rot, "carrier_output": rc}
    for k, ((c0, t0), (c, t)) in enumerate(zip(planet_states(p, 0.0), planet_states(p, psi))):
        move[f"planet_{k}"] = (Pos(c.real, c.imag, 0)
                               * Location((0, 0, 0), (0, 0, math.degrees(t - t0)))
                               * Pos(-c0.real, -c0.imag, 0))  # fmt: skip
        for n in home:
            if n.endswith(f"_{k}") and n.startswith("_bearing"):
                move[n] = Pos(c.real - c0.real, c.imag - c0.imag, 0)
            elif n.endswith(f"_{k}") and n.startswith("_bolt"):
                move[n] = rc
    return {k: move[k] * v if k in move else v for k, v in home.items()}


# large face down; housing_lower on its lip top (the split face)
PRINT_FLIP = {"lid", "carrier_output", "planet_0", "planet_1", "planet_2", "housing_lower"}


def print_pose(name: str, part: Part) -> Part:
    if name in PRINT_FLIP:
        part = part.rotate(Axis.X, 180)
    bb = part.bounding_box()
    return Pos(-bb.center().X, -bb.center().Y, -bb.min.Z) * part


def max_bridge(part: Part) -> float:
    """Widest flat ceiling of a part in print pose: max over downward planar faces above the bed
    of the smaller bounding-box side (a circular ceiling counts as its outer diameter)."""
    out = 0.0
    for f in part.faces():
        if f.geom_type != GeomType.PLANE or f.bounding_box().max.Z < 1e-6:
            continue
        if f.normal_at().Z < -0.999 and f.area > 0.1:
            bb = f.bounding_box()
            out = max(out, min(bb.size.X, bb.size.Y))
    return out


def build(p: Params, out: Path) -> dict[str, Path]:
    """STEP per part + assembly; STL (print pose) per printed part; assembled STLs."""
    out.mkdir(parents=True, exist_ok=True)
    asm_dir = out / "assembled"
    asm_dir.mkdir(exist_ok=True)
    for old in [*out.glob("*.step"), *out.glob("*.stl"), *asm_dir.glob("*.stl")]:
        old.unlink()
    parts = make_parts(p)
    files = {}
    for name, part in parts.items():
        part.label = name.lstrip("_")
        f = out / f"{part.label}.step"
        export_step(part, str(f))
        files[f.name] = f
        if not name.startswith("_"):
            f = out / f"{name}.stl"
            export_stl(print_pose(name, part), str(f))
            files[f.name] = f
        f = asm_dir / f"{part.label}.stl"
        export_stl(part, str(f))
        files[f"assembled/{f.name}"] = f
    export_step(Compound(label="stepped_reducer", children=list(parts.values())),
                str(out / "assembly.step"))  # fmt: skip
    files["assembly.step"] = out / "assembly.step"
    (out / "parts.json").write_text(json.dumps(
        {"printed": [n for n in parts if not n.startswith("_")],
         "assembled": [n.lstrip("_") for n in parts]}, indent=1))  # fmt: skip
    files["parts.json"] = out / "parts.json"
    return files


def bom(p: Params) -> str:
    rows = [
        ("motor iFlight GM5208-24", 1),
        (f"{p.bearing} {BEARINGS[p.bearing]} (in the planets' large gears)", p.n_planets),
    ]
    if p.pin_kind == "m8_bolt":
        rows.append((f"M8x{p.bolt_len} low head (DIN 7984) carrier -> planet pin", p.n_planets))
    rows += [
        ((f"M3x{p.motor_screw_len(p.adapter_t, p.rotor_thread_depth)} ultra-low head "
          f"(ESCO, head {p.head_h}) adapter -> rotor"), p.rotor_n),
        (f"M3x{p.motor_screw_len(p.base_clamp, p.stator_thread_depth)} base -> stator",
         p.stator_n),
        (f"M3x{p.screw_len(p.base_clamp)} base -> housing_lower (hex hole)", p.n_base_screws),
        (f"M3x{p.screw_len(p.joint_clamp)} housing_upper -> housing_lower (hex hole)", p.n_joint),
    ]  # fmt: skip
    return "\n".join(f"  {n:3d} x {s}" for s, n in rows)


def summary(p: Params) -> str:
    z = p.zs_
    s = size_stepped(p)
    lines = [
        (f"teeth sun/large/small/ring={p.zs}/{p.z1}/{p.z2}/{p.zr} m1={p.m1} m2={p.m2} "
         f"x s/1/2/r={p.xs}/{p.x1}/{p.x2}/{p.xr:.3f} ratio={p.ratio:.3f} "
         f"assembly (zs z2 + zr z1)/3={(p.zs * p.z2 + p.zr * p.z1) / p.n_planets:g}"),
        (f"a={p.a:.3f} ring aw={math.degrees(p.aw2):.2f} deg ring tip/root d="
         f"{2 * p.r_ring_tip:.2f}/{2 * p.r_ring_root:.2f} swing R={p.R_sweep:.2f}"),
        (f"body d={s['D_body']:.2f} max d (cap)={s['D_max']:.2f} height={s['height']:.2f} "
         f"(bottom {z['bottom']:.2f} .. out {z['out']:.2f}) b1={p.b1} b2={p.b2} "
         f"hub_d={2 * p.hub_r:.2f} pin={p.pin_kind} d={2 * p.pin_r:.2f} journal bore "
         f"d={2 * p.journal_r:.2f}"),
        "z: " + " ".join(f"{k}={v:.2f}" for k, v in z.items()),
        ("tip land sun/large/small/ring="
         f"{tip_land(p, p.m1, p.zs, p.xs, p.r_tip_sun):.2f}/"
         f"{tip_land(p, p.m1, p.z1, p.x1, p.r_tip1):.2f}/"
         f"{tip_land(p, p.m2, p.z2, p.x2, p.r_tip2):.2f}/{ring_tip_land(p):.2f}"),
    ]  # fmt: skip
    for k in ("R_body from", "R_shell from", "joint_r from", "s0 from", "height terms"):
        lines.append(f"{k}: " + ", ".join(f"{n}={v:.2f}" for n, v in s[k].items()))
    lines.append("slack (>= 0): " + ", ".join(f"{n}={v:.2f}" for n, v in s["slack"].items()))
    for T in (p.torque_cont, p.torque_peak):
        lines.append(f"reference T_in={T} N m: Lewis "
                     + " ".join(f"{g}={v:.2f}" for g, v in lewis_stress(p, T).items())
                     + " MPa; " + ", ".join(f"{g}={v:.2f}"
                                           for g, v in reference_loads(p, T).items()))  # fmt: skip
    lines.append("BOM (bought):\n" + bom(p))
    return "\n".join(lines)


if __name__ == "__main__":
    p = Params()
    print(summary(p))
    for k, v in wall_thicknesses(p).items():
        print(f"  wall {k}: {v:.2f}")
    for w in design_warnings(p) + [f"thin wall {k}: {v:.2f}" for k, v in thin_walls(p).items()]:
        print("WARNING:", w, file=sys.stderr)
    for f in build(p, Path("out_stepped")).values():
        print(f)
