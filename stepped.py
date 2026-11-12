"""Parametric stepped-planet (compound) planetary reducer, ratio 10, for the GM5208-24 (build123d).

Each planet is one print with two gears on a common axis: the large gear meshes with the sun
(herringbone, module m1), the small gear meshes with the fixed ring (spur, module m2).
Ring fixed, sun input, carrier output: ratio = 1 + zr*z1 / (zs*z2).

Axis = Z, z = 0 on the motor rotor's top face. Stack from -Z to +Z:
  base (screwed to the stator from below) | motor in the housing tube | sun adapter (on the rotor)
  | housing lip (floor under the planets) | ring teeth (housing) around the small gears and the
  sun shaft | large gears + sun teeth inside a thin shell | carrier (output) | lid.

Why this layout: the large gears' swing circle (a + large tip + 0.5) is wider than the ring, so
they sit ABOVE the ring and the cluster (sun + planets + carrier) is lowered into the housing
with the small gears entering the one-piece spur ring. No screw passes beside the large gears;
the shell there is min_wall + 0.8 thick. The lid is held by a bayonet (3 lugs in L-slots cut
0.8 mm into the shell) because any axial screw would need a boss inside the swing circle.
The motor (d 60) cannot pass the ring bore (tip d 59.2), so the base is a separate plate and
the housing slides over the motor from above.

Assembly order:
  1 screw sun_adapter to the rotor; 2 screw base to the stator (cable through the floor notch);
  3 lower the housing over motor + adapter, screw it to the base from below (4 M3);
  4 press 685ZZ into both ends of each planet; mesh the 3 planets radially around the sun
    (herringbone cannot slide axially), phase marks outward; put the carrier on top and push
    each Ø5 shoulder bolt up through both bearings into the carrier, screw its M4 end into the
    tapped hole above the socket (dowel alternative: pins pressed into the carrier first);
  5 lower this cluster: the small gears slide into the ring (turn the cluster < 1 ring pitch to
    find the mesh), the sun's hex socket slides onto the adapter boss (turn the rotor by hand);
  6 drop the lid with its lugs in the entry slots, turn it by lock_deg.

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

import planetary as pl
from planetary import _hexagon, half_thickness, herringbone, inv, outline, ring_half_space
from reducer import BEARINGS, MIN, SCREW_LENGTHS, _ring_holes


@dataclass
class Params:
    # gears: sun / planet large (mesh 1, module m1) and planet small / ring (mesh 2, module m2)
    zs: int = 15
    z1: int = 35
    z2: int = 14
    zr: int = 54
    n_planets: int = 3
    m1: float = 0.9
    m2: float = 1.1
    alpha_deg: float = 20.0
    helix_deg: float = 20.0  # sun-large herringbone, as planetary.py
    xs: float = 0.13  # smallest shift against undercut for z = 15
    x1: float = -0.13  # xs + x1 = 0: standard centre distance a = m1 (zs + z1) / 2
    # 0.4 (study: 0.2) keeps 1.2 mm between the small root and the 685 seat; the tip is then
    # trimmed by small_tip_trim to keep the tip land >= tip_min. xr follows from coaxiality.
    x2: float = 0.4
    b1: float = 7.5  # sun-large face width (herringbone total), Lewis <= 10 MPa at 0.4 N m
    b2: float = 13.5  # small-ring face width
    sun_tip_trim: float = 0.1
    small_tip_trim: float = 0.15
    flank_clr: float = 0.05
    fillet_m: float = 0.38
    ring_fillet_m: float = 0.3  # 0.38 makes the two root fillets of the narrow ring space cross
    pts_per_flank: int = 8
    tip_min: float = 0.45
    # strength (Lewis, Kv = 1, ring form factor of a rack 0.154: both conservative)
    torque_cont: float = 0.25
    torque_peak: float = 0.4
    sigma_allow: float = 10.0
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
    # layout
    adapter_t: float = 3.5
    head_h: float = 3.0
    head_gap: float = 0.5
    hex_af: float = 7.0
    hex_engage: float = 5.0
    hex_clr: float = 0.15
    sun_shaft_r: float = 9.0
    gear_gap: float = 0.5  # small gear to lip, ring top to large gear, large gear to lid band
    sweep_clr: float = 0.5  # large tip circle to the shell bore
    shell_margin: float = 0.8  # shell = min_wall + shell_margin
    carrier_gap: float = 1.0
    carrier_t: float = 8.0
    carrier_wall: float = 2.0
    carrier_rim_clr: float = 1.0
    hub_step: float = 2.0
    hub_stub: float = 3.0
    lid_t: float = 3.0
    lid_bore_clr: float = 0.05
    lid_axial_clr: float = 0.05
    lid_rib_h: float = 1.0
    base_t: float = 4.0
    base_clamp: float = 3.5
    wire_notch_deg: float = 45.0
    wire_notch_r0: float = 19.0
    wire_notch_w: float = 10.0
    # planet support: 2 x 685ZZ per planet on a Ø5 shank whose top end sits in a carrier socket.
    # "shoulder_bolt": Ø5 x pin_len shoulder (ISO 7379 style, d1 5 -> M4), head under the small
    # bearing's inner race, M4 end screwed into a hole tapped in the carrier above the socket.
    # "dowel": Ø5 x pin_len steel pin pressed into the same socket (head and thread unused).
    bearing: str = "685ZZ"
    bearing_shoulder: float = 1.0
    pin_kind: str = "shoulder_bolt"
    pin_len: float = 30.0  # shoulder (or dowel) length, s0 .. socket bottom
    pin_head_d: float = 8.0  # socket head of a Ø5/M4 shoulder bolt (vendor tables: 8-9)
    pin_head_h: float = 4.0
    pin_thread_len: float = 5.0  # vendor tables: 4.5-5
    pin_tap_d: float = 3.3  # M4 tap drill; tapped in the print
    pin_socket_d: float = 5.0  # dowel: printed holes come out undersize -> press fit; tune
    pin_socket_min: float = 8.0
    # bayonet lid
    n_lugs: int = 3
    lug_deg: float = 20.0
    lock_deg: float = 30.0
    lug_h: float = 3.0
    lug_depth: float = 0.8  # = shell_margin: leaves min_wall behind the slot
    # fasteners and FDM
    n_base_screws: int = 4
    base_screw_r: float = 36.0
    n_out_holes: int = 4
    m3_clear_d: float = 3.4
    m3_tap_d: float = 2.5
    screw_d: float = 3.0
    screw_engage: float = 6.0
    hex_hole_extra: float = 1.0
    screw_head_d: float = 5.5
    fit_clr: float = 0.05
    slide_clr: float = 0.05
    min_wall: float = 1.2

    # --- derived ---
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
    def R_sweep(self) -> float:  # shell bore around the large gears
        return self.a + self.r_tip1 + self.sweep_clr

    @property
    def R_body(self) -> float:
        return max(self.R_sweep + self.min_wall + self.shell_margin, self.R_motor_bore + 2.0)

    @property
    def R_motor_bore(self) -> float:
        return self.motor_d / 2 + self.motor_clr

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
    def hex_circ(self) -> float:
        return self.screw_d / 2 / math.cos(math.pi / 6)

    @property
    def cbore_r(self) -> float:
        return self.screw_head_d / 2 + self.slide_clr

    @property
    def shoulder_r(self) -> float:
        d, D, _ = BEARINGS[self.bearing]
        return d / 2 + (D - d) / 8

    @property
    def carrier_r(self) -> float:
        return self.a + self.shoulder_r + self.carrier_wall

    @property
    def hub_r(self) -> float:
        return self.carrier_r - self.hub_step

    @property
    def cavity_r(self) -> float:
        return self.carrier_r + self.carrier_rim_clr

    @property
    def zs_(self) -> dict[str, float]:
        """Key heights: s0/s1 small gear (ring teeth ring_lo..s1), L0/L1 large gear."""
        s0 = self.adapter_t + self.head_h + self.head_gap
        s1 = s0 + self.b2
        L0 = s1 + self.gear_gap
        L1 = L0 + self.b1
        c0 = L1 + self.carrier_gap
        c1 = c0 + self.carrier_t
        lid = c1 + self.lid_axial_clr + self.lid_rib_h
        base_top = -self.motor_h
        return {"s0": s0, "ring_lo": s0 - self.gear_gap, "s1": s1, "L0": L0, "L1": L1,
                "Lmid": (L0 + L1) / 2, "band": L1 + self.gear_gap, "c0": c0, "c1": c1,
                "lid": lid, "top": lid + self.lid_t, "out": lid + self.lid_t + self.hub_stub,
                "base_top": base_top, "bottom": base_top - self.base_t}  # fmt: skip

    @property
    def height(self) -> float:
        return self.zs_["out"] - self.zs_["bottom"]

    @property
    def lug_z(self) -> tuple[float, float]:
        """Lug bottom/top; the top bears on the slot ceiling (lid pulled up)."""
        z1 = self.zs_["lid"] - 2.0
        return z1 - self.lug_h, z1

    @property
    def lip_top(self) -> float:  # under the shoulder-bolt heads, gear_gap below them
        return self.zs_["s0"] - self.pin_head_h - self.gear_gap

    @property
    def pin_socket(self) -> float:
        """Pin from the planet bottom (s0) up into the carrier; depth from the boss bottom."""
        z = self.zs_
        return self.pin_len - (z["L1"] - z["s0"])

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
    """Lewis root stress [MPa], Kv = 1. Ring side Ft2 = Ft1 r1/r2 (planet moment balance)."""
    ft1 = torque_in / (p.m1 * p.zs / 2e3) / p.n_planets
    ft2 = ft1 * (p.m1 * p.z1) / (p.m2 * p.z2)

    def y(z: int) -> float:
        return 0.154 - 0.912 / z

    return {"sun": ft1 / (p.b1 * y(p.zs) * math.pi * p.m1),
            "planet_large": ft1 / (p.b1 * y(p.z1) * math.pi * p.m1),
            "planet_small": ft2 / (p.b2 * y(p.z2) * math.pi * p.m2),
            "ring": ft2 / (p.b2 * 0.154 * math.pi * p.m2)}  # fmt: skip


def design_warnings(p: Params) -> list[str]:
    w = []
    if (p.zs * p.z2 + p.zr * p.z1) % p.n_planets:
        w.append("assembly: (zs z2 + zr z1)/N not integer")
    if abs(p.xs + p.x1) > 1e-9:
        w.append("xs + x1 != 0: model assumes the standard sun-large centre distance")
    if not 0.80 < math.cos(p.aw2) < 0.99:
        w.append(f"ring mesh working pressure angle {math.degrees(p.aw2):.1f} deg out of range")
    if planet_tip_gap(p) < p.min_wall:
        w.append(f"planet tips {planet_tip_gap(p):.2f} mm apart")
    for name, m, z, x, ra in (("sun", p.m1, p.zs, p.xs, p.r_tip_sun),
                              ("planet_large", p.m1, p.z1, p.x1, p.r_tip1),
                              ("planet_small", p.m2, p.z2, p.x2, p.r_tip2)):  # fmt: skip
        if x < xmin(p, z) - 1e-9:
            w.append(f"undercut: {name} z={z} x={x}")
        if tip_land(p, m, z, x, ra) < p.tip_min:
            w.append(f"{name} tip land {tip_land(p, m, z, x, ra):.2f} < {p.tip_min}")
    if ring_tip_land(p) < p.tip_min:
        w.append(f"ring tip land {ring_tip_land(p):.2f} < {p.tip_min}")
    if p.r_ring_tip <= p.m2 * p.zr / 2 * math.cos(p.alpha):
        w.append("ring tip inside its base circle")
    w += [f"internal {k} interference" for k, v in internal_interference(p).items() if v]
    for T in (p.torque_cont, p.torque_peak):
        for g, s in lewis_stress(p, T).items():
            if s > p.sigma_allow + 1e-9:
                w.append(f"Lewis: {g} {s:.2f} MPa > {p.sigma_allow} at {T} N m")
    if p.pin_socket < p.pin_socket_min:
        w.append(f"pin socket {p.pin_socket:.1f} < {p.pin_socket_min}")
    if p.zs_["L1"] + p.pin_socket > p.zs_["c1"]:
        w.append("pin socket runs out of the carrier plate")
    if p.a - p.r_tip2 - p.sun_shaft_r < 1.0:
        w.append("sun shaft within 1 mm of the small gear tips")
    if p.lip_top < p.motor_top_clr + p.min_wall:
        w.append("housing lip under the bolt heads thinner than min_wall")
    if p.zs_["L1"] + p.pin_socket + p.pin_thread_len + 1 > p.zs_["out"] - p.min_wall:
        w.append("pin thread hole reaches the output face")
    if p.R_motor_bore < p.r_ring_tip:
        w.append("motor would pass the ring bore: base could be part of the housing")
    if p.lug_z[0] < p.zs_["band"]:
        w.append("bayonet lugs reach into the large-gear zone")
    for face, clamp, depth in (("rotor", p.adapter_t, p.rotor_thread_depth),
                               ("stator", p.base_clamp, p.stator_thread_depth)):  # fmt: skip
        if p.motor_screw_len(clamp, depth) - clamp < 2:
            w.append(f"{face} screws engage < 2 mm (no stock length fits)")
    return w


def wall_thicknesses(p: Params) -> dict[str, float]:
    _, D, B = BEARINGS[p.bearing]
    z = p.zs_
    seat = D / 2 + p.fit_clr
    return {
        "sun shaft around hex socket": p.sun_shaft_r - p.hex_sock,
        "planet small root to 685 seat": p.r_root(p.m2, p.z2, p.x2) - seat,
        "planet large root to 685 seat": p.r_root(p.m1, p.z1, p.x1) - seat,
        "planet web between the 685 seats": (z["L1"] - B) - (z["s0"] + B),
        "ring root to outside": p.R_body - p.r_ring_root,
        "shell around the large-gear swing": p.R_body - p.R_sweep,
        "shell behind the bayonet slot": p.R_body - p.R_sweep - p.lug_depth,
        "housing around motor": p.R_body - p.R_motor_bore,
        "housing lip (axial)": p.lip_top - p.motor_top_clr,
        "housing: base screw hole to motor bore": p.base_screw_r - p.hex_circ - p.R_motor_bore,
        "housing: base screw hole to outside": p.R_body - p.base_screw_r - p.hex_circ,
        "base: counterbore to outside": p.R_body - p.base_screw_r - p.cbore_r,
        "carrier plate around pins": p.carrier_r - p.a - p.pin_socket_d / 2,
        "lid": p.lid_t,
        "lid band around carrier cavity": p.R_sweep - p.slide_clr - p.cavity_r,
        "adapter around screw holes": p.adapter_r - p.rotor_pcd / 2 - p.m3_clear_d / 2,
        "adapter boss around rotor screw heads": p.rotor_pcd / 2 - p.screw_head_d / 2 - p.hex_sock,
        "base under counterbores": p.base_clamp,
        "base: centre hole to stator screw": p.stator_pcd / 2
        - p.m3_clear_d / 2
        - (p.shaft_od / 2 + 0.5),
        "base: wire notch to stator counterbore": p.stator_pcd
        / 2
        * math.sin(
            math.radians(
                min(p.wire_notch_deg % (360 / p.stator_n), -p.wire_notch_deg % (360 / p.stator_n))
            )
        )
        - p.wire_notch_w / 2
        - p.cbore_r,
        "hub around output holes": p.hub_r / 2 - p.m3_tap_d / 2,
        "planet tips": planet_tip_gap(p),
    }


def thin_walls(p: Params) -> dict[str, float]:
    return {k: v for k, v in wall_thicknesses(p).items() if v < p.min_wall - 1e-9}


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


def _sector(r0: float, r1: float, deg0: float, deg1: float, z0: float, h: float) -> Part:
    """Annular sector, angles in degrees."""
    n = 12
    a = np.radians(np.linspace(deg0, deg1, n))
    pts = [(r1 * math.cos(t), r1 * math.sin(t)) for t in a]
    pts += [(r0 * math.cos(t), r0 * math.sin(t)) for t in a[::-1]]
    return _spur(pts, z0, z0 + h)


def lug_angles(p: Params) -> list[float]:
    """Centre angles [deg] of the lugs in the locked position (between the planets)."""
    return [180 / p.n_lugs + 360 * k / p.n_lugs for k in range(p.n_lugs)]


def lid_part(p: Params, turn_deg: float = 0.0) -> Part:
    """Lid; turn_deg = -lock_deg puts the lugs over the entry slots."""
    z = p.zs_
    lid = Pos(0, 0, z["lid"]) * Cylinder(p.R_body, p.lid_t, align=MIN)
    r_band = p.R_sweep - p.slide_clr
    lid += Pos(0, 0, z["band"]) * Cylinder(r_band, z["lid"] - z["band"], align=MIN)
    # carrier cavity under the rib
    lid -= Pos(0, 0, z["band"]) * Cylinder(p.cavity_r, z["lid"] - p.lid_rib_h - z["band"],
                                           align=MIN)  # fmt: skip
    lz0, lz1 = p.lug_z
    for c in lug_angles(p):
        lid += _sector(r_band - 0.5, p.R_sweep + p.lug_depth - p.slide_clr,
                       c - p.lug_deg / 2, c + p.lug_deg / 2, lz0, lz1 - lz0)  # fmt: skip
    g = p.lid_bore_clr
    bore = p.hub_r + g
    lid -= Pos(0, 0, z["c1"]) * Cone(p.carrier_r + g * math.sqrt(2), bore,
                                     p.hub_step + g * (math.sqrt(2) - 1), align=MIN)  # fmt: skip
    lid -= Pos(0, 0, z["c1"]) * Cylinder(bore, z["top"] - z["c1"], align=MIN)
    return lid.rotate(Axis.Z, turn_deg) if turn_deg else lid


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

    adapter = Cylinder(p.adapter_r, p.adapter_t, align=MIN)
    adapter += _hexagon(p.hex_af, p.adapter_t, z["s0"] + p.hex_engage - p.adapter_t)
    adapter -= _ring_holes(p.rotor_n, p.rotor_pcd / 2, p.m3_clear_d, 0, p.adapter_t,
                           math.radians(p.rotor_phase_deg))  # fmt: skip
    parts["sun_adapter"] = adapter

    # sun: herringbone teeth at the large-gear level on a shaft down to the adapter hex
    sun = _herring(p, outline(q1, p.zs, half_thickness(q1, p.zs, p.xs), p.r_tip_sun,
                              p.r_root(p.m1, p.zs, p.xs)), q1.r(p.zs), 1)  # fmt: skip
    zt = z["L0"] - p.gear_gap
    sun += Pos(0, 0, z["s0"]) * Cylinder(p.sun_shaft_r, zt - z["s0"], align=MIN)
    sun += Pos(0, 0, zt) * Cylinder(p.r_root(p.m1, p.zs, p.xs), p.gear_gap, align=MIN)
    parts["sun"] = sun - _hexagon(p.hex_af + 2 * p.hex_clr, z["s0"], p.hex_engage + p.gear_gap)

    # planet: spur small gear s0..L0 + herringbone large gear L0..L1, tooth 0 of both on +X
    small = _spur(outline(q2, p.z2, half_thickness(q2, p.z2, p.x2), p.r_tip2,
                          p.r_root(p.m2, p.z2, p.x2)), z["s0"], z["L0"])  # fmt: skip
    large = _herring(p, outline(q1, p.z1, half_thickness(q1, p.z1, p.x1), p.r_tip1,
                                p.r_root(p.m1, p.z1, p.x1)), q1.r(p.z1), -1)  # fmt: skip
    planet = small.fuse(large).clean()
    planet -= Pos(0, 0, z["L1"] - B) * Cylinder(D / 2 + p.fit_clr, B, align=MIN)
    planet -= Pos(0, 0, z["s0"]) * Cylinder(D / 2 + p.fit_clr, B, align=MIN)
    planet -= Pos(0, 0, z["s0"]) * Cylinder(D / 2 - p.bearing_shoulder, z["L1"] - z["s0"],
                                            align=MIN)  # fmt: skip
    mark_r = (D / 2 + p.r_root(p.m1, p.z1, p.x1)) / 2  # phase mark on the large gear's top face
    planet -= Pos(mark_r, 0, z["L1"] - 0.4) * Cylinder(0.6, 0.4, align=MIN)

    carrier = Pos(0, 0, z["c0"]) * Cylinder(p.carrier_r, p.carrier_t, align=MIN)
    carrier += Pos(0, 0, z["c1"]) * Cone(p.carrier_r, p.hub_r, p.hub_step, align=MIN)
    zh = z["c1"] + p.hub_step
    carrier += Pos(0, 0, zh) * Cylinder(p.hub_r, z["out"] - zh, align=MIN)
    carrier -= _ring_holes(p.n_out_holes, p.hub_r / 2, p.m3_tap_d, z["out"] - 8, 8)
    pin = f"_pin_{p.pin_kind}_5x{p.pin_len:g}"
    shank = Cylinder(d / 2, p.pin_len, align=MIN)
    if p.pin_kind == "shoulder_bolt":
        shank += Pos(0, 0, -p.pin_head_h) * Cylinder(p.pin_head_d / 2, p.pin_head_h, align=MIN)
        shank += Pos(0, 0, p.pin_len) * Cylinder(
            p.pin_tap_d / 2, p.pin_thread_len, align=MIN
        )  # M4 at its core d
    for k, (c, turn) in enumerate(planet_states(p, 0.0)):
        x, y = c.real, c.imag
        parts[f"planet_{k}"] = Pos(x, y, 0) * planet.rotate(Axis.Z, math.degrees(turn))
        for tag, zb in (("large", z["L1"] - B), ("small", z["s0"])):
            parts[f"_bearing_{p.bearing}_{tag}_{k}"] = Pos(x, y, zb) * (
                Cylinder(D / 2, B, align=MIN) - Cylinder(d / 2, B, align=MIN)
            )
        parts[f"{pin}_{k}"] = Pos(x, y, z["s0"]) * shank
        carrier += Pos(x, y, z["L1"]) * Cylinder(p.shoulder_r, p.carrier_gap, align=MIN)
        carrier -= Pos(x, y, z["L1"]) * Cylinder(p.pin_socket_d / 2, p.pin_socket, align=MIN)
        carrier -= Pos(x, y, z["L1"] + p.pin_socket) * Cylinder(
            p.pin_tap_d / 2, p.pin_thread_len + 1, align=MIN
        )
    parts["carrier_output"] = carrier

    # base: under the stator, screwed to it and to the housing from below
    base = Pos(0, 0, z["bottom"]) * Cylinder(p.R_body, p.base_t, align=MIN)
    base -= Pos(0, 0, z["bottom"]) * Cylinder(p.shaft_od / 2 + 0.5, p.base_t, align=MIN)
    sph = math.radians(p.stator_phase_deg)
    cb = p.base_t - p.base_clamp
    base -= _ring_holes(p.stator_n, p.stator_pcd / 2, p.m3_clear_d, z["bottom"], p.base_t, sph)
    base -= _ring_holes(p.stator_n, p.stator_pcd / 2, 2 * p.cbore_r, z["bottom"], cb, sph)
    base -= _ring_holes(p.n_base_screws, p.base_screw_r, p.m3_clear_d, z["bottom"], p.base_t, sph)
    base -= _ring_holes(p.n_base_screws, p.base_screw_r, 2 * p.cbore_r, z["bottom"], cb, sph)
    base -= (Rot(0, 0, p.stator_phase_deg + p.wire_notch_deg)
             * Pos(p.wire_notch_r0, -p.wire_notch_w / 2, z["bottom"])
             * Box(p.R_body - p.wire_notch_r0 + 1, p.wire_notch_w, p.base_t, align=MIN))  # fmt: skip
    parts["base"] = base

    # housing: motor tube, lip, spur ring, thin shell with bayonet slots
    h = Pos(0, 0, z["base_top"]) * Cylinder(p.R_body, z["lid"] - z["base_top"], align=MIN)
    h -= Pos(0, 0, z["base_top"]) * Cylinder(p.R_motor_bore, p.motor_h + p.motor_top_clr,
                                             align=MIN)  # fmt: skip
    h -= Cylinder(p.lip_r, z["s1"], align=MIN)
    h -= Pos(0, 0, p.lip_top) * Cylinder(p.r_ring_tip, z["ring_lo"] - p.lip_top, align=MIN)
    qr = replace(q2, fillet_m=p.ring_fillet_m)
    ring = _spur(outline(qr, p.zr, ring_half_space(qr), p.r_ring_root, p.r_ring_tip, cutter=True),
                 z["ring_lo"], z["s1"])  # fmt: skip
    h -= ring.rotate(Axis.Z, math.degrees(ring_phase(p)))
    h -= Pos(0, 0, z["s1"]) * Cylinder(p.R_sweep, z["lid"] - z["s1"], align=MIN)
    lz0, lz1 = p.lug_z
    r_slot = p.R_sweep + p.lug_depth
    for c in lug_angles(p):
        w = p.lug_deg / 2 + 1.0  # 1 deg angular clearance
        h -= _sector(p.R_sweep - 1, r_slot, c - p.lock_deg - w, c - p.lock_deg + w, lz0 - 0.5,
                     z["lid"] - lz0 + 0.5)  # entry slot, open at the top  # fmt: skip
        h -= _sector(p.R_sweep - 1, r_slot, c - p.lock_deg - w, c + w, lz0 - 0.5,
                     lz1 + p.lid_axial_clr - lz0 + 0.5)  # locking pocket  # fmt: skip
    hd = p.screw_len(p.base_clamp) - p.base_clamp + p.hex_hole_extra
    for k in range(p.n_base_screws):
        ph = sph + 2 * math.pi * k / p.n_base_screws
        h -= Pos(p.base_screw_r * math.cos(ph), p.base_screw_r * math.sin(ph), 0) * _hexagon(
            p.screw_d, z["base_top"], hd
        )
    parts["housing"] = h
    parts["lid"] = lid_part(p)
    return parts


def pose(p: Params, home: dict[str, Part], psi: float) -> dict[str, Part]:
    """Sun at angle psi [rad] (ring fixed, carrier psi/ratio)."""
    rot = Location((0, 0, 0), (0, 0, math.degrees(psi)))
    move = {"sun": rot, "sun_adapter": rot,
            "carrier_output": Location((0, 0, 0), (0, 0, math.degrees(psi / p.ratio)))}  # fmt: skip
    for k, ((c0, t0), (c, t)) in enumerate(zip(planet_states(p, 0.0), planet_states(p, psi))):
        move[f"planet_{k}"] = (Pos(c.real, c.imag, 0)
                               * Location((0, 0, 0), (0, 0, math.degrees(t - t0)))
                               * Pos(-c0.real, -c0.imag, 0))  # fmt: skip
        shift = Pos(c.real - c0.real, c.imag - c0.imag, 0)
        for n in home:
            if n.endswith(f"_{k}") and n.startswith(("_bearing", "_pin")):
                move[n] = shift
    return {k: move[k] * v if k in move else v for k, v in home.items()}


PRINT_FLIP = {"lid", "carrier_output"}


def print_pose(name: str, part: Part) -> Part:
    if name in PRINT_FLIP:
        part = part.rotate(Axis.X, 180)
    bb = part.bounding_box()
    return Pos(-bb.center().X, -bb.center().Y, -bb.min.Z) * part


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
        (f"{p.bearing} {BEARINGS[p.bearing]}", 2 * p.n_planets),
        (
            f"shoulder bolt d5 x {p.pin_len:g} M4 x {p.pin_thread_len:g} (or dowel d5 x {p.pin_len:g})",
            p.n_planets,
        ),
        (f"M3x{p.motor_screw_len(p.adapter_t, p.rotor_thread_depth)} adapter -> rotor", p.rotor_n),
        (f"M3x{p.motor_screw_len(p.base_clamp, p.stator_thread_depth)} base -> stator", p.stator_n),
        (f"M3x{p.screw_len(p.base_clamp)} base -> housing (hex hole)", p.n_base_screws),
    ]
    return "\n".join(f"  {n:3d} x {s}" for s, n in rows)


def summary(p: Params) -> str:
    z = p.zs_
    lines = [
        (
            f"teeth sun/large/small/ring={p.zs}/{p.z1}/{p.z2}/{p.zr} m1={p.m1} m2={p.m2} "
            f"x s/1/2/r={p.xs}/{p.x1}/{p.x2}/{p.xr:.3f} ratio={p.ratio:.3f} "
            f"assembly (zs z2 + zr z1)/3={(p.zs * p.z2 + p.zr * p.z1) / p.n_planets:g}"
        ),
        (
            f"a={p.a:.3f} ring aw={math.degrees(p.aw2):.2f} deg ring tip/root d="
            f"{2 * p.r_ring_tip:.2f}/{2 * p.r_ring_root:.2f} swing R={p.R_sweep:.2f}"
        ),
        (
            f"outer_d={2 * p.R_body:.2f} height={p.height:.2f} (bottom {z['bottom']:.2f} .. out "
            f"{z['out']:.2f}) b1={p.b1} b2={p.b2} hub_d={2 * p.hub_r:.2f} pin socket={p.pin_socket}"
        ),
        "z: " + " ".join(f"{k}={v:.2f}" for k, v in z.items()),
        (
            "tip land sun/large/small/ring="
            f"{tip_land(p, p.m1, p.zs, p.xs, p.r_tip_sun):.2f}/"
            f"{tip_land(p, p.m1, p.z1, p.x1, p.r_tip1):.2f}/"
            f"{tip_land(p, p.m2, p.z2, p.x2, p.r_tip2):.2f}/{ring_tip_land(p):.2f}"
        ),
    ]
    for T in (p.torque_cont, p.torque_peak):
        lines.append(f"Lewis T_in={T} N m: "
                     + " ".join(f"{g}={v:.2f}" for g, v in lewis_stress(p, T).items()))  # fmt: skip
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
