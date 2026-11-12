"""Parametric 1-stage planetary reducer (herringbone) for the iFlight GM5208-24, built with build123d.

Axis = Z, z = 0 on the motor rotor's top face. Stack from -Z to +Z:
  housing base (screwed to the stator from below) | motor | sun adapter (screwed to the rotor)
  | sun + 3 planets (608ZZ inside each) between the lower ring half (part of the housing) and
  the upper ring half (part of ring_lid) | carrier (output) | lid with the output hub bore.
Ring fixed, sun input, carrier output: ratio = 1 + Zr/Zs (Singh2025_COMPAct p2 eq. 2).

Design basis (literature/surveys/2026-11-10_planetary-qdd-design-basis.md; pages are PDF pages):
  concentric Zr = Zs + 2Zp, assembly (Zs+Zr)/N integer, adjacency Zp+2 < (Zs+Zp) sin(180/N)
      -- Khk2025_InternalRingGearsCatalog p3 (Application Hints 3)
  internal gear interference limits (pinion teeth vs ring teeth) -- KHK catalog p2 table
  undercut-free shift x >= 1 - z sin^2(alpha)/2 -- textbook formula (survey 2.2, not checked
      against KHK because that PDF could not be downloaded)
  KHK Gear Technical Reference "Calculation of Gear Dimensions" (web, no pages;
      https://khkgears.net/new/gear_knowledge/gear_technical_reference/calculation_gear_dimensions.html):
      Table 4.3 shifted spur pair (x1+x2 = 0 -> alpha_w = alpha, a = m(z1+z2)/2; z1=12 needs a
      positive shift against undercut), Table 4.6 internal gear ha2 = (1-x2)m, da2 = d2 - 2 ha2,
      interference: involute (standard 20 deg ring z2 > 34), trochoid (standard z2 - z1 > 9),
      Tables 4.9/4.11 helical: we use transverse module/pressure angle (db = d cos alpha_t)
  modified Lewis sigma = Ft/(b y Kv P) -- Singh2025_ExternalInternalPlanetary p3 eq. 9;
      y = 0.154 - 0.912/z and Kv = 6.1/(6.1+V) are the usual textbook forms (the paper takes
      them from its ref. [22], which we did not have)
  helix angle 20 deg -- not stated in any source; derived from Nguyen2025 p4 Table 2: helical/
      herringbone gears had 47 teeth vs 50 spur teeth (m=1) at the same 50 mm centre distance,
      cos(beta) = 47/50 -> beta = 19.9 deg (the other modules give 18-23 deg)

Assembly order (why the sun floats and the ring is split at the herringbone apex):
  1 motor into housing, screw stator from below; 2 screw sun_adapter to rotor;
  3 put the 3 planets (bearings pressed in) on the carrier pins, put the sun between them
    (herringbone gears can be joined radially, not axially);
  4 lower this cluster into the housing: the planets' lower halves screw into the lower ring
    half while the cluster turns as one body; the sun's hex socket slides onto the adapter boss
    (turn the rotor by hand to align);
  5 lower ring_lid: its ring half screws onto the planets' upper halves; screw it down with M3.
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
    Location,
    Part,
    Polyline,
    Pos,
    RegularPolygon,
    Rot,
    Solid,
    Wire,
    export_step,
    export_stl,
    extrude,
)

from reducer import BEARINGS, MIN, SCREW_LENGTHS, _ring_holes

# KHK catalog p2: ring teeth -> (min pinion teeth: involute interference,
# max: trochoid interference, max: trimming interference). Valid for KHK standard (x=0) gears.
KHK_INTERNAL = {60: (21, 51, 43), 80: (20, 72, 64), 100: (19, 92, 84),
                120: (19, 112, 104), 160: (19, 152, 144), 200: (18, 192, 184)}  # fmt: skip
BUILD_VOLUME = (260.0, 260.0, 300.0)  # Creality Hi


@dataclass
class Params:
    # gears (transverse module; teeth: sun, planet, ring; N planets)
    zs: int = 12
    zp: int = 24
    zr: int = 60
    n_planets: int = 3
    m: float = 1.2  # smallest module whose Zp=24 planet still holds a 608ZZ with min_wall
    alpha_deg: float = 20.0
    helix_deg: float = 20.0  # Nguyen2025 p4 Table 2 (derived, see module docstring)
    face_w: float = 12.0  # total herringbone width (two 6 mm halves)
    xs: float = 0.3  # sun shift: z=12 needs x >= 0.30 against undercut
    xp: float = -0.3  # planet shift; xs + xp = 0 keeps the standard centre distance; xr = xp
    # backlash: each flank moved back along its normal (per flank), and planet centre distance
    flank_clr: float = 0.05  # same as the ball reducer's 0.05 mm per side (assumption)
    center_clr: float = 0.0  # +: more backlash at the sun, LESS at the ring (keep 0)
    # strength check
    torque_cont: float = 0.25  # N m, GM5208-24 load torque (iFlight page: 0.18-0.25)
    torque_peak: float = 0.4  # N m, short-time assumption (max torque not published)
    # MPa, placeholder to be set by break tests: gear-test limits are 38-72 % of VDI 2736 values
    # even for moulded PA6 (Pogacnik2015 p10 Table 7); stronger catalogue materials broke first
    # in printed planetaries (Yoshida2023 p4-5)
    sigma_allow: float = 10.0
    input_rpm: float = 1000.0  # assumption, only enters Kv
    # motor: iFlight product page + RobotShop GM5208-24 drawing
    # https://cdn.robotshop.com/media/i/ifl/rb-ifl-11/img/gm5208-24.webp
    # rotor top: 12-M3 on PCD 25 (4), 30 (4), outer 4 at 45 deg; stator bottom: 15-M3 on PCD 25,
    # 30, 44. Hollow shaft 20.5 long < body 23.2 -> no protrusion assumed. Phase angles read
    # from the drawing (holes on the X/Y axes) -> check on the part.
    motor_d: float = 60.0  # envelope over the measured 59.3 (iFlight page: 59.5)
    motor_h: float = 23.25  # measured stator face to rotor top 23.2-23.3 (iFlight page: 31.2)
    shaft_od: float = 15.0
    shaft_id: float = 12.6
    # hollow shaft turns with the rotor (checked by the user); protrusion not measured yet
    shaft_protrude_top: float = 0.0
    shaft_protrude_bottom: float = 0.0
    rotor_n: int = 4
    rotor_pcd: float = 25.0  # alternative: 30.0
    rotor_phase_deg: float = 0.0
    stator_n: int = 4
    stator_pcd: float = 30.0  # alternatives: 25.0, 44.0 (7 holes, one missing at the connector)
    stator_phase_deg: float = 0.0
    # Section view of the RobotShop drawing, pixel estimate (4.70 px/mm from 23.2, 4.65 from
    # 59.5; +-0.45 mm): rotor top = left face (bell end plate), plate and through holes
    # 15 px = 3.2 mm; stator bottom = right face, full-diameter base, blind holes 15 px = 3.2 mm.
    # Usable depth = 3.2 - 0.5 margin. Neither face has a protrusion.
    rotor_thread_depth: float = 2.7
    stator_thread_depth: float = 2.7
    motor_clr: float = 0.75  # radial, rotating bell to housing bore
    # layout
    adapter_t: float = 3.5  # M3x6 then engages 2.5 mm <= rotor_thread_depth
    head_h: float = 3.0  # M3 socket head height
    head_gap: float = 0.5
    hex_af: float = 7.0  # sun adapter boss across flats
    hex_engage: float = 5.0
    gear_gap: float = 0.5  # planet faces to the ring steps above and below
    carrier_gap: float = 1.0  # planet top to carrier plate (pin shoulder height)
    carrier_t: float = 4.0
    carrier_wall: float = 2.0  # carrier plate beyond the pin shoulders
    carrier_rim_clr: float = 1.0
    hub_step: float = 2.0  # carrier plate radius - hub radius, bridged by a 45 deg cone
    hub_stub: float = 3.0
    lid_t: float = 3.0
    lid_bore_clr: float = 0.05
    lid_axial_clr: float = 0.05
    lid_rib_h: float = 1.0  # lid underside below the lid plate (over the carrier cavity)
    base_t: float = 4.0
    # 0: the drawing shows a full-diameter stator base under the bell (bell ends ~5 mm above it)
    base_boss_h: float = 0.0
    base_clamp: float = 3.5  # base under the stator screw heads (rest is counterbore)
    # motor cable leaves the stator bottom straight down at ~24 mm, midway between PCD-30 holes
    wire_notch_deg: float = 45.0  # relative to stator_phase_deg
    wire_notch_r0: float = 19.0  # radial notch through the floor, from r0 out through the wall
    wire_notch_w: float = 10.0
    bearing: str = "608ZZ"
    bearing_lip_t: float = 1.5
    bearing_shoulder: float = 1.0
    # fasteners and FDM
    n_ring_screws: int = 4
    n_out_holes: int = 4
    m3_clear_d: float = 3.4
    m3_tap_d: float = 2.5
    screw_d: float = 3.0  # inscribed diameter of the hex hole the screw cuts into
    screw_engage: float = 6.0
    hex_hole_extra: float = 1.0
    screw_head_d: float = 5.5
    fit_clr: float = 0.05
    slide_clr: float = 0.05
    min_wall: float = 1.2
    nozzle_d: float = 0.4
    pts_per_flank: int = 8
    # root fillet radius / m: 0.38 = tip radius of the standard basic rack (customary value,
    # not found on the KHK pages we read); against root cracks seen in Yoshida2023 p4-5
    fillet_m: float = 0.38
    sun_tip_trim: float = 0.1  # mm off the sun tip radius so its tip land stays >= nozzle_d

    # --- derived ---
    @property
    def alpha(self) -> float:
        return math.radians(self.alpha_deg)

    @property
    def xr(self) -> float:
        return self.xp

    @property
    def ratio(self) -> float:
        return 1 + self.zr / self.zs

    @property
    def a(self) -> float:  # sun-planet centre distance
        return self.m * (self.zs + self.zp) / 2 + self.center_clr

    def r(self, z: int) -> float:
        return self.m * z / 2

    @property
    def r_ring_tip(self) -> float:  # KHK internal: d_a = d - 2m(1 - x)
        return self.r(self.zr) - self.m * (1 - self.xr)

    @property
    def r_ring_root(self) -> float:
        return self.r(self.zr) + self.m * (1.25 + self.xr)

    def r_tip(self, z: int, x: float) -> float:
        return self.r(z) + self.m * (1 + x)

    @property
    def r_tip_sun(self) -> float:
        return self.r_tip(self.zs, self.xs) - self.sun_tip_trim

    def r_root(self, z: int, x: float) -> float:
        return self.r(z) - self.m * (1.25 - x)

    @property
    def hex_circ(self) -> float:
        return self.screw_d / 2 / math.cos(math.pi / 6)

    @property
    def cbore_r(self) -> float:
        return self.screw_head_d / 2 + self.slide_clr

    @property
    def R_s(self) -> float:  # ring screw circle: 1.1 * min_wall outside the ring root
        return self.r_ring_root + 1.1 * self.min_wall + max(self.m3_clear_d / 2, self.hex_circ)

    @property
    def R_body(self) -> float:
        return self.R_s + self.cbore_r + 1.1 * self.min_wall

    @property
    def R_motor_bore(self) -> float:
        return self.motor_d / 2 + self.motor_clr

    @property
    def pin_r(self) -> float:
        return BEARINGS[self.bearing][0] / 2 - self.fit_clr

    @property
    def shoulder_r(self) -> float:  # stays on the inner race
        d, D, _ = BEARINGS[self.bearing]
        return d / 2 + (D - d) / 8

    @property
    def carrier_r(self) -> float:
        return self.a + self.shoulder_r + self.carrier_wall

    @property
    def hub_r(self) -> float:
        return self.carrier_r - self.hub_step

    @property
    def adapter_r(self) -> float:
        return self.rotor_pcd / 2 + self.screw_head_d / 2 + 1.1 * self.min_wall

    @property
    def zs_(self) -> dict[str, float]:
        """Key heights."""
        g0 = self.adapter_t + self.head_h + self.head_gap  # gear bottom
        gt = g0 + self.face_w
        c0 = gt + self.carrier_gap
        lid = c0 + self.carrier_t + self.lid_axial_clr + self.lid_rib_h
        base_top = -self.motor_h - self.base_boss_h
        return {"g0": g0, "mid": g0 + self.face_w / 2, "gt": gt, "c0": c0, "c1": c0 + self.carrier_t,
                "lid": lid, "top": lid + self.lid_t, "out": lid + self.lid_t + self.hub_stub,
                "base_top": base_top, "bottom": base_top - self.base_t}  # fmt: skip

    @property
    def height(self) -> float:
        return self.zs_["out"] - self.zs_["bottom"]

    @property
    def ring_clamp(self) -> float:  # ring_lid under the screw heads (ring teeth zone)
        return self.face_w / 2 + self.gear_gap

    def screw_len(self, clamp: float) -> int:
        return next(L for L in SCREW_LENGTHS if L >= clamp + self.screw_engage)

    def motor_screw_len(self, clamp: float, depth: float) -> int:
        """Longest stock screw that does not go deeper than depth into the motor (0: none fits)."""
        return max((L for L in SCREW_LENGTHS if L <= clamp + depth), default=0)


# ---------------------------------------------------------------- gear geometry
def inv(a: float) -> float:
    return math.tan(a) - a


def half_thickness(p: Params, z: int, x: float) -> float:
    """Half tooth angle at the pitch circle [rad], flank_clr removed (external gear)."""
    s = p.m * (math.pi / 2 + 2 * x * math.tan(p.alpha))
    rb = p.r(z) * math.cos(p.alpha)
    return s / (2 * p.r(z)) - p.flank_clr / rb


def ring_half_space(p: Params) -> float:
    """Half space angle of the ring at the pitch circle [rad] (KHK: s2 = m(pi/2 - 2x tan a))."""
    e = p.m * (math.pi / 2 + 2 * p.xr * math.tan(p.alpha))
    return e / (2 * p.r(p.zr)) + p.flank_clr / (p.r(p.zr) * math.cos(p.alpha))


def _fillet(pts: list, corners: set[int], rho: float, n: int = 4) -> list:
    """Replace each corner point by an arc of radius rho tangent to its two neighbours' lines.

    Points closer to the corner than the tangent length are dropped (flank points near the root).
    """
    P = np.array(pts)
    drop, arcs = set(), {}
    for i in corners:
        c = P[i]
        u, v = P[i - 1] - c, P[(i + 1) % len(P)] - c
        u, v = u / np.linalg.norm(u), v / np.linalg.norm(v)
        th = math.acos(float(np.clip(u @ v, -1, 1)))
        t = rho / math.tan(th / 2)
        for step in (-1, 1):  # drop neighbours inside the tangent length
            k = i + step
            while np.linalg.norm(P[k % len(P)] - c) < t and (k % len(P)) not in corners:
                drop.add(k % len(P))
                k += step
        a, b = c + u * t, c + v * t
        m = c + (u + v) / np.linalg.norm(u + v) * rho / math.sin(th / 2)
        a0, a1 = math.atan2(*(a - m)[::-1]), math.atan2(*(b - m)[::-1])
        a1 = a0 + (a1 - a0 + math.pi) % (2 * math.pi) - math.pi  # short way round
        arcs[i] = [tuple(m + rho * np.array([math.cos(q), math.sin(q)]))
                   for q in np.linspace(a0, a1, n + 1)]  # fmt: skip
    out = []
    for i, q in enumerate(pts):
        if i in arcs:
            out += arcs[i]
        elif i not in drop:
            out.append(q)
    return out


def outline(
    p: Params, z: int, hp: float, r_tip: float, r_root: float, cutter: bool = False
) -> list[tuple[float, float]]:
    """Closed tooth outline (tooth j centred on angle 2 pi j/z). hp = half angle at pitch.

    For the ring this is the cutter: its 'teeth' are the ring's tooth spaces, tip = ring root.
    Below the base circle the flank continues radially.
    """
    rb = p.r(z) * math.cos(p.alpha)

    def h(R: float) -> float:
        return hp + inv(p.alpha) - inv(math.acos(min(1.0, rb / R)))

    fl = [(R, h(R)) for R in np.linspace(max(rb, r_root), r_tip, p.pts_per_flank)]
    if r_root < rb:
        fl.insert(0, (r_root, h(rb)))
    pts = []
    for j in range(z):
        c = 2 * math.pi * j / z
        pts += [(R * math.cos(c - a), R * math.sin(c - a)) for R, a in fl]
        pts += [(R * math.cos(c + a), R * math.sin(c + a)) for R, a in reversed(fl)]
    # fillet the gear's root corners; for the ring cutter its tip corners are the ring's root
    n = len(fl)
    local = (n - 1, n) if cutter else (0, 2 * n - 1)
    return _fillet(pts, {2 * n * j + k for j in range(z) for k in local}, p.fillet_m * p.m)


def tip_land(p: Params, z: int, x: float, ra: float | None = None) -> float:
    """Chord width of the tooth tip [mm] at tip radius ra (default r_tip)."""
    rb, ra = p.r(z) * math.cos(p.alpha), ra or p.r_tip(z, x)
    return 2 * ra * (half_thickness(p, z, x) + inv(p.alpha) - inv(math.acos(rb / ra)))


def herringbone(p: Params, pts, r_pitch: float, hand: int, z_lo: float, z_hi: float) -> Solid:
    """Twisted extrusion with the apex at the gear mid-plane; hand=+1 sun, -1 planet and ring.

    The section at height z is turned by hand*tan(helix)/r_pitch*|z - z_mid| (rad), so slices of
    sun, planets and ring at the same z are a valid carrier-fixed mesh (see test_planetary).
    """
    zm = p.zs_["mid"]
    rate = hand * math.degrees(math.tan(math.radians(p.helix_deg)) / r_pitch)  # deg per mm
    w = Wire(Polyline(*pts, close=True).edges())
    out = []
    if z_lo < zm:  # lower half: turn rate*(zm - z_lo) at z_lo, 0 at zm
        sec = Pos(0, 0, z_lo) * w.rotate(Axis.Z, rate * (zm - z_lo))
        out.append(Solid.extrude_linear_with_rotation(sec, (0, 0, z_lo), (0, 0, zm - z_lo),
                                                      -rate * (zm - z_lo)))  # fmt: skip
    if z_hi > zm:
        sec = Pos(0, 0, zm) * w
        out.append(Solid.extrude_linear_with_rotation(sec, (0, 0, zm), (0, 0, z_hi - zm),
                                                      rate * (z_hi - zm)))  # fmt: skip
    return out[0] if len(out) == 1 else out[0].fuse(out[1]).clean()


# ---------------------------------------------------------------- kinematics
def planet_phase(p: Params, k: int) -> tuple[float, float]:
    """(carrier angle of planet k, its own turn [rad]) at input 0, so its teeth fill sun spaces."""
    th = 2 * math.pi * k / p.n_planets
    u = (p.zs * th / (2 * math.pi)) % 1  # sun phase at th (0 = tooth centre)
    return th, th + math.pi - 2 * math.pi * (0.5 - u) / p.zp


def ring_phase(p: Params) -> float:
    """Turn [rad] of the ring cutter (whose teeth are the ring's spaces) so that the ring teeth
    fill planet 0's spaces; the half pitch moves from a ring tooth centre to a space centre."""
    _, phi = planet_phase(p, 0)
    v = p.zp * (0 - phi) / (2 * math.pi)
    return -2 * math.pi * (v - 0.5) / p.zr + math.pi / p.zr


def planet_states(p: Params, psi: float) -> list[tuple[complex, float]]:
    """(centre, absolute turn [rad]) of each planet when the sun has turned psi (ring fixed)."""
    pc = psi / p.ratio
    out = []
    for k in range(p.n_planets):
        th, phi = planet_phase(p, k)
        out.append((p.a * np.exp(1j * (th + pc)), phi + pc - (psi - pc) * p.zs / p.zp))
    return out


# ---------------------------------------------------------------- checks
def khk_internal_limits(zr: int) -> tuple[int, int, int] | None:
    # shortcut: nearest table row below zr, upper limits shifted by zr - row; interpolate if
    # a ring between rows becomes the design
    rows = [k for k in KHK_INTERNAL if k <= zr]
    if not rows:
        return None
    k = max(rows)
    lo, troch, trim = KHK_INTERNAL[k]
    return lo, troch + zr - k, trim + zr - k


def design_warnings(p: Params) -> list[str]:
    w = []
    if p.zr != p.zs + 2 * p.zp:
        w.append(f"concentric: Zr={p.zr} != Zs + 2Zp = {p.zs + 2 * p.zp}")
    if (p.zs + p.zr) % p.n_planets:
        w.append(f"assembly: (Zs+Zr)/N = {(p.zs + p.zr) / p.n_planets} not integer")
    if not p.zp + 2 < (p.zs + p.zp) * math.sin(math.pi / p.n_planets):
        w.append("adjacency: Zp+2 >= (Zs+Zp) sin(180/N)")
    if planet_tip_gap(p) < p.min_wall:
        w.append(f"planet tips {planet_tip_gap(p):.2f} mm apart")
    for name, z, x, ra in (("sun", p.zs, p.xs, p.r_tip_sun), ("planet", p.zp, p.xp, None)):
        xmin = 1 - z * math.sin(p.alpha) ** 2 / 2
        if x < xmin - 1e-9:
            w.append(f"undercut: {name} z={z} x={x} < {xmin:.3f}")
        if tip_land(p, z, x, ra) < p.nozzle_d:
            w.append(f"{name} tip land {tip_land(p, z, x, ra):.2f} < nozzle {p.nozzle_d}")
    if abs(p.xs + p.xp) > 1e-9:
        w.append("xs + xp != 0: model assumes the standard centre distance")
    lim = khk_internal_limits(p.zr)
    if lim is None:
        w.append(f"KHK internal table starts at 60 teeth (Zr={p.zr})")
    else:
        lo, troch, trim = lim
        if not lo <= p.zp <= min(troch, trim):
            w.append(f"KHK internal: Zp={p.zp} outside {lo}..{min(troch, trim)}")
    if p.zr <= 34 or p.zr - p.zp <= 9:  # KHK dimension page, Table 4.6 notes (standard gears)
        w.append(f"KHK internal: needs Zr > 34 and Zr - Zp > 9 (Zr={p.zr}, Zp={p.zp})")
    if p.r_ring_tip <= p.r(p.zr) * math.cos(p.alpha):
        w.append("involute interference: ring tip inside its base circle")
    for T in (p.torque_cont, p.torque_peak):
        for gear, s in lewis_stress(p, T).items():
            if s > p.sigma_allow:
                w.append(f"Lewis: {gear} {s:.1f} MPa > {p.sigma_allow} at {T} N m")
    if p.shaft_protrude_top + p.fit_clr >= p.adapter_t:
        w.append("shaft protrudes through the sun adapter")
    for face, clamp, depth in (("rotor", p.adapter_t, p.rotor_thread_depth),
                               ("stator", p.base_clamp, p.stator_thread_depth)):  # fmt: skip
        if p.motor_screw_len(clamp, depth) - clamp < 2:
            w.append(f"{face} screws engage < 2 mm (no stock length fits)")
    return w


def planet_tip_gap(p: Params) -> float:
    return 2 * p.a * math.sin(math.pi / p.n_planets) - 2 * p.r_tip(p.zp, p.xp)


def lewis_stress(p: Params, torque_in: float) -> dict[str, float]:
    """Modified Lewis root stress [MPa] (Singh2025_ExternalInternalPlanetary p3 eq. 9).

    Ft from equal load sharing over N planets at the sun pitch circle; same Ft on the planet.
    Transverse values and the full herringbone width (no helix credit).
    """
    rs = p.r(p.zs)
    ft = torque_in / (rs / 1000) / p.n_planets  # N
    v = p.input_rpm * 2 * math.pi / 60 * rs / 1000  # m/s at the sun pitch circle
    kv = 6.1 / (6.1 + v)
    P = math.pi * p.m
    return {g: ft / (p.face_w * (0.154 - 0.912 / z) * kv * P)
            for g, z in (("sun", p.zs), ("planet", p.zp))}  # fmt: skip


def wall_thicknesses(p: Params) -> dict[str, float]:
    D = BEARINGS[p.bearing][1]
    hex_sock = (p.hex_af / 2 + p.fit_clr) / math.cos(math.pi / 6)
    return {
        "sun root to hex socket": p.r_root(p.zs, p.xs) - hex_sock,
        "planet root to bearing seat": p.r_root(p.zp, p.xp) - (D / 2 + p.fit_clr),
        "planet lip under bearing": p.face_w - BEARINGS[p.bearing][2],
        "ring root to screw hole": p.R_s - max(p.m3_clear_d / 2, p.hex_circ) - p.r_ring_root,
        "body outside counterbore": p.R_body - p.R_s - p.cbore_r,
        "housing around motor": p.R_body - p.R_motor_bore,
        "carrier plate around pins": p.carrier_r - p.a - p.pin_r,
        "lid": p.lid_t,
        "adapter around screw holes": p.adapter_r - p.rotor_pcd / 2 - p.m3_clear_d / 2,
        "adapter boss around rotor screw heads": p.rotor_pcd / 2 - p.screw_head_d / 2 - hex_sock,
        "base under counterbores": p.base_clamp,
        "base: centre hole to stator screw": p.stator_pcd / 2
        - p.m3_clear_d / 2
        - (p.shaft_od / 2 + 0.5),
        # perpendicular distance from the notch centre line (lower bound on the true gap)
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
        **(
            {
                "adapter: shaft pocket to rotor screw": p.rotor_pcd / 2
                - p.m3_clear_d / 2
                - (p.shaft_od / 2 + p.fit_clr)
            }
            if p.shaft_protrude_top > 0
            else {}
        ),
        "planet tips": planet_tip_gap(p),
    }


def thin_walls(p: Params) -> dict[str, float]:
    return {k: v for k, v in wall_thicknesses(p).items() if v < p.min_wall}


def candidates() -> list[dict]:
    """Ratio candidates at a common 72 mm ring pitch diameter (comparable outer size)."""
    out = []
    for zs, zp, zr in ((12, 24, 60), (12, 30, 72), (12, 36, 84)):
        q = replace(Params(), zs=zs, zp=zp, zr=zr, m=round(72 / zr, 3))
        out.append({"teeth": (zs, zp, zr), "m": q.m, "ratio": q.ratio, "outer_d": 2 * q.R_body,
                    "sigma_peak_sun": lewis_stress(q, q.torque_peak)["sun"],
                    "608_wall": wall_thicknesses(q)["planet root to bearing seat"],
                    "reflected_inertia_x": q.ratio**2})  # fmt: skip
    return out


# ---------------------------------------------------------------- solids
def _hexagon(af: float, z0: float, h: float) -> Part:
    return Pos(0, 0, z0) * extrude(RegularPolygon(af / 2, 6, major_radius=False), h)


def _bearing(p: Params, z0: float) -> Part:
    d, D, B = BEARINGS[p.bearing]
    return Pos(0, 0, z0) * (Cylinder(D / 2, B, align=MIN) - Cylinder(d / 2, B, align=MIN))


def ring_cutter(p: Params, z_lo: float, z_hi: float) -> Solid:
    pts = outline(p, p.zr, ring_half_space(p), p.r_ring_root, p.r_ring_tip, cutter=True)
    return herringbone(p, pts, p.r(p.zr), -1, z_lo, z_hi).rotate(
        Axis.Z, math.degrees(ring_phase(p))
    )


def make_parts(p: Params) -> dict[str, Part]:
    """Parts at input angle 0, planets at their home angles. Names starting with '_' are bought."""
    z = p.zs_
    _, D, B = BEARINGS[p.bearing]
    parts: dict[str, Part] = {}

    motor = Pos(0, 0, -p.motor_h) * Cylinder(p.motor_d / 2, p.motor_h, align=MIN)
    z_sh = -p.motor_h - p.shaft_protrude_bottom
    motor += Pos(0, 0, z_sh) * Cylinder(p.shaft_od / 2, z_sh * -1 + p.shaft_protrude_top, align=MIN)
    parts["_motor_GM5208-24"] = motor - Pos(0, 0, z_sh) * Cylinder(
        p.shaft_id / 2, -z_sh + p.shaft_protrude_top, align=MIN
    )

    adapter = Cylinder(p.adapter_r, p.adapter_t, align=MIN)
    adapter += _hexagon(p.hex_af, p.adapter_t, z["g0"] + p.hex_engage - p.adapter_t)
    adapter -= _ring_holes(
        p.rotor_n, p.rotor_pcd / 2, p.m3_clear_d, 0, p.adapter_t, math.radians(p.rotor_phase_deg)
    )
    if p.shaft_protrude_top > 0:  # shaft turns with the adapter: a close pocket is enough
        adapter -= Cylinder(p.shaft_od / 2 + p.fit_clr, p.shaft_protrude_top + p.fit_clr, align=MIN)
    parts["sun_adapter"] = adapter

    sun = herringbone(
        p, outline(p, p.zs, half_thickness(p, p.zs, p.xs), p.r_tip_sun,
                   p.r_root(p.zs, p.xs)), p.r(p.zs), 1, z["g0"], z["gt"],
    )  # fmt: skip
    parts["sun"] = sun - _hexagon(p.hex_af + 2 * p.fit_clr, z["g0"], p.hex_engage + p.gear_gap)

    planet = herringbone(
        p, outline(p, p.zp, half_thickness(p, p.zp, p.xp), p.r_tip(p.zp, p.xp),
                   p.r_root(p.zp, p.xp)), p.r(p.zp), -1, z["g0"], z["gt"],
    )  # fmt: skip
    planet -= Pos(0, 0, z["gt"] - B) * Cylinder(D / 2 + p.fit_clr, B, align=MIN)
    planet -= Pos(0, 0, z["g0"]) * Cylinder(D / 2 - p.bearing_shoulder, p.face_w, align=MIN)
    carrier = Pos(0, 0, z["c0"]) * Cylinder(p.carrier_r, p.carrier_t, align=MIN)
    # 45 deg cone from the plate to the hub: prints face down without support
    carrier += Pos(0, 0, z["c1"]) * Cone(p.carrier_r, p.hub_r, p.hub_step, align=MIN)
    zh = z["c1"] + p.hub_step
    carrier += Pos(0, 0, zh) * Cylinder(p.hub_r, z["out"] - zh, align=MIN)
    carrier -= _ring_holes(p.n_out_holes, p.hub_r / 2, p.m3_tap_d, z["out"] - 8, 8)
    for k, (c, turn) in enumerate(planet_states(p, 0.0)):
        x, y = c.real, c.imag
        parts[f"planet_{k}"] = Pos(x, y, 0) * planet.rotate(Axis.Z, math.degrees(turn))
        parts[f"_bearing_{p.bearing}_{k}"] = Pos(x, y, 0) * _bearing(p, z["gt"] - B)
        carrier += Pos(x, y, z["gt"] - B) * Cylinder(p.pin_r, B, align=MIN)
        carrier += Pos(x, y, z["gt"]) * Cylinder(p.shoulder_r, p.carrier_gap, align=MIN)
    parts["carrier_output"] = carrier

    # housing: base under the stator, tube around the motor, lower ring half on top
    ring_lo = z["g0"] - p.gear_gap
    h = Pos(0, 0, z["bottom"]) * Cylinder(p.R_body, z["mid"] - z["bottom"], align=MIN)
    h -= Pos(0, 0, z["base_top"]) * Cylinder(p.R_motor_bore, ring_lo - z["base_top"], align=MIN)
    boss_r = p.stator_pcd / 2 + p.screw_head_d / 2 + 1.1 * p.min_wall
    if p.base_boss_h > 0:
        h += Pos(0, 0, z["base_top"]) * Cylinder(boss_r, p.base_boss_h, align=MIN)
    h -= ring_cutter(p, ring_lo, z["mid"])
    # through hole for the turning hollow shaft (0.5 mm radial clearance)
    h -= Pos(0, 0, z["bottom"]) * Cylinder(p.shaft_od / 2 + 0.5, p.base_t + p.base_boss_h,
                                           align=MIN)  # fmt: skip
    stator_t = p.base_t + p.base_boss_h
    sph = math.radians(p.stator_phase_deg)
    h -= _ring_holes(p.stator_n, p.stator_pcd / 2, p.m3_clear_d, z["bottom"], stator_t, sph)
    h -= _ring_holes(p.stator_n, p.stator_pcd / 2, 2 * p.cbore_r, z["bottom"],
                     stator_t - p.base_clamp, sph)  # fmt: skip
    hd = p.screw_len(p.ring_clamp) - p.ring_clamp + p.hex_hole_extra
    for k in range(p.n_ring_screws):
        ph = 2 * math.pi * (k + 0.5) / p.n_ring_screws
        h -= Pos(p.R_s * math.cos(ph), p.R_s * math.sin(ph), 0) * _hexagon(
            p.screw_d, z["mid"] - hd, hd
        )
    # motor cable exits the stator bottom straight down: radial notch through the whole floor
    nl = p.R_body - p.wire_notch_r0 + 1
    h -= (
        Rot(0, 0, p.stator_phase_deg + p.wire_notch_deg)
        * Pos(p.wire_notch_r0, -p.wire_notch_w / 2, z["bottom"])
        * Box(nl, p.wire_notch_w, stator_t, align=MIN)
    )
    parts["housing"] = h

    # ring_lid: upper ring half, cavity for the carrier plate, lid with the hub bore and rib
    z_ledge = z["gt"] + p.gear_gap
    rl = Pos(0, 0, z["mid"]) * Cylinder(p.R_body, z["top"] - z["mid"], align=MIN)
    rl -= ring_cutter(p, z["mid"], z_ledge)
    rl -= Pos(0, 0, z_ledge) * Cylinder(p.carrier_r + p.carrier_rim_clr, z["lid"] - z_ledge,
                                        align=MIN)  # fmt: skip
    # lid underside spans the cavity; bore + 45 deg cone parallel to the hub, lid_bore_clr normal
    rl += Pos(0, 0, z["lid"] - p.lid_rib_h) * Cylinder(p.carrier_r + p.carrier_rim_clr,
                                                      p.lid_rib_h, align=MIN)  # fmt: skip
    g = p.lid_bore_clr
    bore = p.hub_r + g
    rl -= Pos(0, 0, z["c1"]) * Cone(p.carrier_r + g * math.sqrt(2), bore,
                                    p.hub_step + g * (math.sqrt(2) - 1), align=MIN)  # fmt: skip
    rl -= Pos(0, 0, z["c1"]) * Cylinder(bore, z["top"] - z["c1"], align=MIN)
    ph0 = math.pi / p.n_ring_screws
    rl -= _ring_holes(p.n_ring_screws, p.R_s, p.m3_clear_d, z["mid"], z["top"] - z["mid"], ph0)
    zc = z["mid"] + p.ring_clamp
    rl -= _ring_holes(p.n_ring_screws, p.R_s, 2 * p.cbore_r, zc, z["top"] - zc, ph0)
    parts["ring_lid"] = rl
    return parts


def pose(p: Params, home: dict[str, Part], psi: float) -> dict[str, Part]:
    """Move the parts to sun angle psi [rad] (ring fixed, carrier psi/ratio)."""
    move = {"sun": Location((0, 0, 0), (0, 0, math.degrees(psi))),
            "sun_adapter": Location((0, 0, 0), (0, 0, math.degrees(psi))),
            "carrier_output": Location((0, 0, 0), (0, 0, math.degrees(psi / p.ratio)))}  # fmt: skip
    for k, ((c0, t0), (c, t)) in enumerate(zip(planet_states(p, 0.0), planet_states(p, psi))):
        # planet k: undo its home placement, turn, move to the new centre
        loc = (Pos(c.real, c.imag, 0) * Location((0, 0, 0), (0, 0, math.degrees(t - t0)))
               * Pos(-c0.real, -c0.imag, 0))  # fmt: skip
        move[f"planet_{k}"] = loc
        move[f"_bearing_{p.bearing}_{k}"] = Pos(c.real - c0.real, c.imag - c0.imag, 0)
    return {k: move[k] * v if k in move else v for k, v in home.items()}


PRINT_FLIP = {"ring_lid", "carrier_output"}  # turned over to avoid overhangs


def print_pose(name: str, part: Part) -> Part:
    if name in PRINT_FLIP:
        part = part.rotate(Axis.X, 180)
    bb = part.bounding_box()
    return Pos(-bb.center().X, -bb.center().Y, -bb.min.Z) * part


def build(p: Params, out: Path) -> dict[str, Path]:
    """STEP for every part + assembly; STL (print pose) for printed parts; assembled STLs."""
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
    export_step(Compound(label="planetary_reducer", children=list(parts.values())),
                str(out / "assembly.step"))  # fmt: skip
    files["assembly.step"] = out / "assembly.step"
    (out / "parts.json").write_text(json.dumps(
        {"printed": [n for n in parts if not n.startswith("_")],
         "assembled": [n.lstrip("_") for n in parts]}, indent=1))  # fmt: skip
    files["parts.json"] = out / "parts.json"
    return files


def summary(p: Params) -> str:
    z = p.zs_
    lines = [
        (
            f"teeth Zs/Zp/Zr={p.zs}/{p.zp}/{p.zr} N={p.n_planets} m={p.m} alpha={p.alpha_deg} "
            f"helix={p.helix_deg} xs/xp/xr={p.xs}/{p.xp}/{p.xr} ratio={p.ratio:.3f}"
        ),
        (
            f"pitch d sun/planet/ring={2 * p.r(p.zs):.2f}/{2 * p.r(p.zp):.2f}/{2 * p.r(p.zr):.2f} "
            f"a={p.a:.3f} ring tip/root d={2 * p.r_ring_tip:.2f}/{2 * p.r_ring_root:.2f}"
        ),
        (
            f"outer_d={2 * p.R_body:.2f} height={p.height:.2f} (base {z['bottom']:.2f} .. output face "
            f"{z['out']:.2f}; motor top z=0) face_w={p.face_w} carrier_r={p.carrier_r:.2f} "
            f"hub_d={2 * p.hub_r:.2f}"
        ),
        (
            f"planet tip gap={planet_tip_gap(p):.2f} tip land sun/planet="
            f"{tip_land(p, p.zs, p.xs, p.r_tip_sun):.2f}/{tip_land(p, p.zp, p.xp):.2f}"
        ),
        (
            f"screws: ring M3x{p.screw_len(p.ring_clamp)} x{p.n_ring_screws}, rotor "
            f"M3x{p.motor_screw_len(p.adapter_t, p.rotor_thread_depth)} x{p.rotor_n}, stator "
            f"M3x{p.motor_screw_len(p.base_clamp, p.stator_thread_depth)} x{p.stator_n}"
        ),
    ]
    for T in (p.torque_cont, p.torque_peak):
        s = lewis_stress(p, T)
        lines.append(f"Lewis T_in={T} N m (T_out ideal {T * p.ratio:.2f}): "
                     + " ".join(f"{g}={v:.2f} MPa ({v / p.sigma_allow:.2f} of {p.sigma_allow})"
                                for g, v in s.items()))  # fmt: skip
    for c in candidates():
        lines.append(
            "candidate "
            + " ".join(f"{k}={v:.2f}" if isinstance(v, float) else f"{k}={v}" for k, v in c.items())
        )
    return "\n".join(lines)


if __name__ == "__main__":
    p = Params()
    print(summary(p))
    for k, v in wall_thicknesses(p).items():
        print(f"  wall {k}: {v:.2f}")
    for w in design_warnings(p) + [f"thin wall {k}: {v:.2f}" for k, v in thin_walls(p).items()]:
        print("WARNING:", w, file=sys.stderr)
    for f in build(p, Path("out_planetary")).values():
        print(f)
