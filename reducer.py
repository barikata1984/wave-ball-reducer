"""Parametric cycloid ball reducer (Wang & Ping 2019, IRJAES 4(2) 185-187), built with build123d.

Axis = Z. Stack from -Z to +Z: 3 end cover | 4 balls | 5 centre disc | 4 balls |
6 hypocycloid disc | 7 balls + 9 cage | 8 epicycloid disc (output).
1 eccentric input shaft runs through 3 and 5 and carries 2 bearing inside 6.
10 housing tube + 11 lid; 3 and 11 are each screwed into hex holes in 10 with short M3 screws + nut traps.
FDM printing: balls and bearing are bought (STEP only); printed parts get STL in print pose.

Conventions taken from the paper:
  Rb = nb*r0                              (Table II row 7)
  K  = O0M/r0 = A/r0                      (eq. 5; Table II's K = A*nb/(2Rb) is read as a typo)
  K2 = Rb*sin(pi/nb)/rb in [1.0, 1.15]    (Table II row 11)
  epicycloid groove (disc 8, ne = nb-1 waves):  z = Rb e^{it} - A e^{i nb t}      (eq. 7/9)
  hypocycloid groove (disc 6, nh = nb+1 waves): z = Rb e^{it} + A e^{-i nb t}     (eq. 8)
  groove depth h = rb - rb*sin(beta)      (eq. 10), cross-section = double circle (gothic arch)

Reduction ratio (derived here, see reduction_ratio()/crossing_counts()): with 6 translating
(eccentric angle psi) and 8 rotating (alpha), the radial gap between the grooves is to first
order A[cos((nb+1)t) + cos((nb-1)(t-alpha)) + cos(t-psi)], which factorises into
A cos(t-psi)[1 + 2cos(nb t + psi)] only when alpha = -2 psi/(nb-1). Hence i = -(nb-1)/2.
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
from build123d import (
    Align,
    Axis,
    Circle,
    Compound,
    Cylinder,
    Location,
    Part,
    Plane,
    Pos,
    RegularPolygon,
    Sphere,
    Spline,
    Vector,
    export_step,
    export_stl,
    extrude,
    sweep,
)

TABLE_I = [
    (11, 0.42, 0.50),
    (28, 0.42, 0.63),
    (49, 0.52, 0.75),
    (70, 0.68, 0.79),
    (90, 0.71, 0.80),
]
# Deep-groove ball bearings, ZZ widths: name -> (bore d, outer D, width B) [mm]
BEARINGS = {
    "683ZZ": (3, 7, 3), "684ZZ": (4, 9, 4), "685ZZ": (5, 11, 5), "686ZZ": (6, 13, 5),
    "687ZZ": (7, 14, 5), "688ZZ": (8, 16, 5), "689ZZ": (9, 17, 5),
    "6800ZZ": (10, 19, 5), "6801ZZ": (12, 21, 5), "6803ZZ": (17, 26, 5),
    "6804ZZ": (20, 32, 7), "6805ZZ": (25, 37, 7), "6806ZZ": (30, 42, 7), "608ZZ": (8, 22, 7),
}  # fmt: skip
MIN = (None, None, Align.MIN)


SCREW_LENGTHS = (6, 8, 12, 16)  # stock M3 socket-head lengths [mm]


@dataclass
class Params:
    rb: float = 3.0  # steel ball radius (both sections), bought part
    nb: int = 21  # reduction balls -> ratio (nb-1)/2 = 10
    # K, K2: smallest housing with A >= 0.6 mm and 10 % curvature margin, see search_K_K2()
    K: float = 0.295  # cycloid factor A/r0 (Table I wants 0.42-0.63, which undercuts)
    K2: float = 2.13  # ball-spacing coefficient, sets Rb (paper: 1.0-1.15)
    beta_deg: float = 32.0  # eq.10 depth parameter; sets disc gap 2 rb sin(beta) -> cage thickness
    gothic_deg: float = 45.0  # contact angle of the gothic arch arcs
    conformity: float = 1.04  # gothic arc radius / groove radius
    groove_floor: float = 3.0  # material under the deepest groove (>= min_wall)
    bearing_lip_t: float = 1.5  # axial lip over a bearing seat
    bearing_boss_wall: float = 2.0  # radial wall of the cover's bearing boss
    rim: float = 4.0  # material outside the outermost groove
    bearing: str = "608ZZ"  # 2: eccentric bearing in disc 6, also 13 input bearing
    # disc 8 hub runs in the lid bore (plain bearing); a rib under the lid takes its lift-off
    hub_r: float = 20.0
    lid_bore_clr: float = 0.05  # radial, hub to lid bore
    lid_axial_clr: float = 0.05  # disc 8 top face to the lid rib
    lid_rib_w: float = 1.5
    lid_rib_h: float = 1.0
    # FDM tolerances (radial, i.e. per side)
    ball_clr: float = 0.05  # ball vs groove / cage pocket
    fit_clr: float = 0.05  # bearing seats and journals in bearing bores
    slide_clr: float = 0.05  # sliding / rotating printed surfaces
    min_wall: float = 1.2
    # design targets
    a_min: float = 0.6  # minimum eccentricity (FDM error 0.1-0.2 mm)
    curv_margin: float = 1.10  # min centre-line curvature radius / groove half-width at the face
    cage_min_t: float = 2.5
    # fasteners and drive
    hex_af: float = 3.0  # hex socket in the input shaft (3 mm hex key / hex bit)
    hex_depth: float = 8.0
    shaft_stub: float = 10.0  # input shaft length below the end cover
    housing_wall: float = 5.0  # body wall; screws sit in outer ears
    lid_t: float = 4.0
    hub_stub: float = 3.0  # output hub protrusion above the lid
    n_screws: int = 4
    m3_clear_d: float = 3.4
    m3_tap_d: float = 2.5  # self-tapping holes in the output hub for attaching the load
    screw_d: float = 3.0  # inscribed diameter of the hex hole in 10 the screw cuts into
    screw_engage: float = 6.0  # minimum thread engagement in the hex hole
    hex_hole_extra: float = 1.0  # hex hole depth beyond the screw tip
    screw_head_d: float = 5.5  # M3 socket head
    clamp_t: float = 4.0  # cover under the head if a counterbore is needed
    n_out_holes: int = 4
    bearing_shoulder: float = 1.0  # radial lip under bearings
    cage_ring: float = 2.0  # cage material inside/outside the ball pockets
    n_oldham: int = 4  # balls per Oldham layer
    pts_per_wave: int = 24  # spline resolution of groove centre lines

    # --- derived ---
    @property
    def Rb(self) -> float:
        return self.K2 * self.rb / math.sin(math.pi / self.nb)

    @property
    def r0(self) -> float:
        return self.Rb / self.nb

    @property
    def A(self) -> float:  # eccentricity, eq. 5
        return self.K * self.r0

    @property
    def h(self) -> float:  # groove depth, eq. 10
        return self.rb - self.rb * math.sin(math.radians(self.beta_deg))

    @property
    def s(self) -> float:  # ball centre height above a grooved face
        return self.rb - self.h

    @property
    def rg(self) -> float:  # groove section radius incl. clearance
        return self.rb + self.ball_clr

    @property
    def groove_depth(
        self,
    ) -> float:  # deepest cut below a grooved face (gothic arch bottom)
        arc = self.conformity * self.rg
        off, g = arc - self.rg, math.radians(self.gothic_deg)
        return math.sqrt(arc**2 - (off * math.cos(g)) ** 2) - off * math.sin(g) - self.s

    @property
    def face_half_width(self) -> float:
        """Groove half-width at the disc face (widest point of the cut; the lens is widest above
        the face). If the centre-line curvature radius is below this, the sweep overlaps itself
        on the inner side of the bend and cuts away the land next to the groove."""
        arc = self.conformity * self.rg
        off, g = arc - self.rg, math.radians(self.gothic_deg)
        return math.sqrt(arc**2 - (self.s + off * math.sin(g)) ** 2) - off * math.cos(g)

    @property
    def contact_offset(self) -> float:
        """Horizontal distance of a gothic contact point from the centre line (rg cos g)."""
        return self.rg * math.cos(math.radians(self.gothic_deg))

    @property
    def cage_t(self) -> float:
        return 2 * self.s - 2 * self.slide_clr

    @property
    def R_out(self) -> float:
        return self.Rb + self.A + self.rb + self.rim

    @property
    def z_top(self) -> float:  # lid bottom = housing top
        return self.zs["8"] + self.t8 + self.lid_axial_clr + self.lid_rib_h

    # plate thicknesses: groove depth on each grooved face + groove_floor. Bearings 2 and 13 are
    # taller; they hang out of disc 6 into disc 5's centre hole and out of the cover in a boss.
    @property
    def t3(self) -> float:  # cover: Oldham slots on top
        return self.groove_depth + self.groove_floor

    @property
    def t5(self) -> float:  # centre disc: slots on both faces
        return 2 * self.groove_depth + self.groove_floor

    @property
    def t6(self) -> float:  # slots below, cycloid groove on top
        return 2 * self.groove_depth + self.groove_floor

    @property
    def t8(self) -> float:  # cycloid groove below; the hub on top is extra
        return self.groove_depth + self.groove_floor

    @property
    def zb13(self) -> float:  # input bearing bottom, in a boss under the cover
        return self.t3 - self.bearing_lip_t - BEARINGS[self.bearing][2]

    @property
    def z_bottom(self) -> float:  # input shaft end
        return min(self.zb13, 0.0) - self.shaft_stub

    @property
    def height(self) -> float:  # shaft end to output face
        return self.z_top + self.lid_t + self.hub_stub - self.z_bottom

    @property
    def zs(self) -> dict[str, float]:
        """Bottom faces of 3, 5, 6, 8 and ball-centre heights 4b, 4a, 7."""
        s = self.s
        z5 = self.t3 + 2 * s
        z6 = z5 + self.t5 + 2 * s
        z8 = z6 + self.t6 + 2 * s
        return {"3": 0.0, "4b": self.t3 + s, "5": z5, "4a": z5 + self.t5 + s, "6": z6,
                "7": z6 + self.t6 + s, "8": z8}  # fmt: skip

    @property
    def cover_clamp(self) -> float:  # cover under the head; counterbored only if no screw fits
        fits = self.t3 + self.screw_engage <= max(SCREW_LENGTHS)
        return self.t3 if fits else self.clamp_t

    def screw_len(self, clamp: float) -> int:
        """Shortest available length with screw_engage beyond a clamped part of this thickness."""
        return next(L for L in SCREW_LENGTHS if L >= clamp + self.screw_engage)

    def hex_depth_10(self, clamp: float) -> float:  # hex hole depth in the housing
        return self.screw_len(clamp) - clamp + self.hex_hole_extra

    @property
    def hex_circ(self) -> float:  # circumradius of the hex hole
        return self.screw_d / 2 / math.cos(math.pi / 6)

    @property
    def cbore_r(self) -> float:
        return self.screw_head_d / 2 + self.slide_clr

    @property
    def R_h(self) -> float:  # housing bore
        return self.R_out + self.A + self.slide_clr

    @property
    def R_cover(self) -> float:  # body radius of cover, housing and lid (ears excluded)
        return self.R_h + self.housing_wall

    @property
    def R_s(self) -> float:  # screw circle: 1.1 * min_wall left between hole and housing bore
        return self.R_h + 1.1 * self.min_wall + max(self.m3_clear_d / 2, self.hex_circ)

    @property
    def ear_r(self) -> float:  # screw ear radius: 1.1 * min_wall around the counterbore
        return self.cbore_r + 1.1 * self.min_wall

    @property
    def cam_r(self) -> float:
        return BEARINGS[self.bearing][0] / 2 - self.fit_clr

    @property
    def shaft_r(self) -> float:  # concentric shaft below the cam, journal in bearing 13
        return self.cam_r

    @property
    def ratio(self) -> float:
        return reduction_ratio(self.nb)


def reduction_ratio(nb: int) -> float:
    """Input turns per output turn (6 non-rotating, 8 output); sign = output opposite to input."""
    return -(nb - 1) / 2


# ---------------------------------------------------------------- curves & checks
def groove_xy(p: Params, kind: str, t: np.ndarray) -> np.ndarray:
    """Complex groove centre line in the disc's own frame. kind: 'epi' (8) or 'hypo' (6)."""
    e = p.Rb * np.exp(1j * t)
    return e - p.A * np.exp(1j * p.nb * t) if kind == "epi" else e + p.A * np.exp(-1j * p.nb * t)


def min_curvature_radius(p: Params, kind: str, n: int = 20000) -> float:
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    dt = t[1]
    z = groove_xy(p, kind, t)
    d1 = (np.roll(z, -1) - np.roll(z, 1)) / (2 * dt)
    d2 = (np.roll(z, -1) - 2 * z + np.roll(z, 1)) / dt**2
    kappa = np.abs(np.imag(np.conj(d1) * d2)) / np.abs(d1) ** 3
    return float(1 / kappa.max())


def design_warnings(p: Params) -> list[str]:
    w = []
    for kind in ("epi", "hypo"):
        rho = min_curvature_radius(p, kind)
        if rho <= p.curv_margin * p.face_half_width:
            w.append(
                f"undercut: {kind} min curvature radius {rho:.3f} "
                f"< {p.curv_margin} * face half-width {p.face_half_width:.3f}"
            )
    if p.A < p.a_min:
        w.append(f"eccentricity A={p.A:.3f} < {p.a_min}")
    if p.rg * math.sin(math.radians(p.gothic_deg)) <= p.s:
        w.append("gothic contact point lies above the disc face")
    if not 1.0 <= p.K2 <= 1.15:
        w.append(f"K2={p.K2:.3f} outside 1.0-1.15")
    lo, hi = next(((lo, hi) for n, lo, hi in TABLE_I if p.nb <= n), (None, None))
    if lo is not None and not lo <= p.K <= hi:
        w.append(f"K={p.K} outside Table I range {lo}-{hi} for nb={p.nb}")
    return w


def _polar_radius(z: np.ndarray, grid: np.ndarray) -> np.ndarray:
    ang = np.mod(np.angle(z), 2 * np.pi)
    o = np.argsort(ang)
    return np.interp(grid, ang[o], np.abs(z)[o], period=2 * np.pi)


def groove_gap(p: Params, psi: float, alpha: float, grid: np.ndarray) -> np.ndarray:
    """Radial gap (6 minus 8) between groove centre lines in the fixed frame vs polar angle."""
    t = np.linspace(0, 2 * np.pi, 40 * p.nb * 8)
    z8 = np.exp(1j * alpha) * groove_xy(p, "epi", t)
    z6 = groove_xy(p, "hypo", t) + p.A * np.exp(1j * psi)
    return _polar_radius(z6, grid) - _polar_radius(z8, grid)


def crossing_counts(p: Params, ratio: float, steps: int = 97) -> list[int]:
    """Number of groove crossings over one input turn when alpha = psi/ratio."""
    grid = np.linspace(0, 2 * np.pi, 8000, endpoint=False)
    out = []
    for psi in np.linspace(0, 2 * np.pi, steps):
        s = np.sign(groove_gap(p, psi, psi / ratio, grid))
        out.append(int(np.sum(s != np.roll(s, 1))))
    return out


def search_K_K2(p: Params, step: float = 0.005) -> tuple[float, float, float]:
    """Smallest outer diameter with A >= a_min and curvature >= curv_margin*face_half_width.

    At fixed K the groove shape scales with K2 (Rb, r0, A all proportional to K2), so for each
    K the minimum K2 is the largest of the two linear bounds. Returns (K, K2, outer diameter).
    """
    best = None
    for K in np.arange(0.05, 0.6 + step / 2, step):
        unit = replace(p, K=float(K), K2=1.0)
        rho1 = min(min_curvature_radius(unit, k, 4000) for k in ("epi", "hypo"))
        K2 = max(1.0, p.a_min / unit.A, p.curv_margin * p.face_half_width / rho1)
        K2 = math.ceil(K2 * 100) / 100
        d = 2 * replace(p, K=float(K), K2=K2).R_cover
        if best is None or d < best[2]:
            best = (round(float(K), 3), K2, d)
    assert best is not None
    return best


def _d(p: Params, kind: str, t: float) -> complex:
    e = 1j * p.Rb * np.exp(1j * t)
    return (
        e - 1j * p.nb * p.A * np.exp(1j * p.nb * t)
        if kind == "epi"
        else e - 1j * p.nb * p.A * np.exp(-1j * p.nb * t)
    )


def ball7_positions(p: Params, psi: float = 0.0) -> np.ndarray:
    """Groove crossings (ball 7 centres, fixed frame) at input angle psi, alpha = psi/ratio.

    Start from the first-order roots nb*t + psi = 2pi/3 + 2pi k and solve
    e^{i alpha} E(t8) = H(t6) + A e^{i psi} by Newton. The geometry does not depend on the
    turning direction: the same crossing is followed whether psi grows or shrinks.
    """
    alpha = psi / p.ratio
    rot, off = np.exp(1j * alpha), p.A * np.exp(1j * psi)
    out = []
    for k in range(p.nb):
        th = (2 * np.pi / 3 + 2 * np.pi * k - psi) / p.nb
        t8, t6 = th - alpha, th
        for _ in range(30):
            F = rot * groove_xy(p, "epi", np.array(t8)) - groove_xy(p, "hypo", np.array(t6)) - off
            J1, J2 = rot * _d(p, "epi", t8), -_d(p, "hypo", t6)
            M = np.array([[J1.real, J2.real], [J1.imag, J2.imag]])
            dt = np.linalg.solve(M, [-F.real, -F.imag])
            t8, t6 = t8 + dt[0], t6 + dt[1]
            if abs(dt).max() < 1e-13:
                break
        out.append(complex(rot * groove_xy(p, "epi", np.array(t8))))
    return np.array(out)


# ---------------------------------------------------------------- printability
def wall_thicknesses(p: Params) -> dict[str, float]:
    """Thinnest printed walls [mm] (closed-form, conservative)."""
    _, D, _ = BEARINGS[p.bearing]
    hex_circ = (p.hex_af / 2 + p.fit_clr) / math.cos(math.pi / 6)
    return {
        "cage web between ball pockets": 2 * p.Rb * math.sin(math.pi / p.nb)
        - 2 * cage_pockets(p, 12)[1],
        "cage ring beside pockets": p.cage_ring,
        "shaft around hex socket": p.shaft_r - hex_circ,
        "disc 6 between bearing seat and groove": p.Rb - p.A - p.rg - (D / 2 + p.fit_clr),
        "groove floor (3, 5, 6, 8)": p.groove_floor,
        "lip over bearings (3, 6)": p.bearing_lip_t,
        "cover bearing boss": p.bearing_boss_wall,
        "housing body wall": p.housing_wall,
        "housing between hex hole and bore": p.R_s - p.hex_circ - p.R_h,
        "ear around hex hole": p.ear_r - p.hex_circ,
        "ear around counterbore": p.ear_r - p.cbore_r,
        "housing between the two hex holes": p.z_top
        - p.t3
        - p.hex_depth_10(p.cover_clamp)
        - p.hex_depth_10(p.lid_t),
    }


def thin_walls(p: Params) -> dict[str, float]:
    return {k: v for k, v in wall_thicknesses(p).items() if v < p.min_wall}


# ---------------------------------------------------------------- solids
def _disc(z0: float, r: float, t: float, hole_r: float = 0.0) -> Part:
    d = Pos(0, 0, z0) * Cylinder(r, t, align=MIN)
    if hole_r:
        d -= Pos(0, 0, z0) * Cylinder(hole_r, t, align=MIN)
    return d


def _gothic_profile(p: Params, plane: Plane, down: int):
    """Double-circle section. down=+1: material below the ball centre in plane-y."""
    arc = p.conformity * p.rg
    g = math.radians(p.gothic_deg)
    off = arc - p.rg
    c1 = plane * Pos(-off * math.cos(g), down * off * math.sin(g)) * Circle(arc)
    c2 = plane * Pos(off * math.cos(g), down * off * math.sin(g)) * Circle(arc)
    return (c1 & c2).face()


def _cycloid_groove(p: Params, kind: str, zc: float, down: int) -> Part:
    n = p.pts_per_wave * (p.nb + 1)
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    z = groove_xy(p, kind, t)
    path = Spline(*[(c.real, c.imag, zc) for c in z], periodic=True)
    p0, t0 = path @ 0, path % 0
    prof = _gothic_profile(p, Plane(origin=p0, x_dir=Vector(0, 0, 1).cross(t0), z_dir=t0), down)
    return sweep(prof, path)


def _slot(p: Params, x: float, y: float, zc: float, along_x: bool) -> Part:
    """Straight ball groove (capsule) allowing +-A travel, centred on a ball."""
    L = 2 * p.A + 2 * p.slide_clr
    body = Cylinder(p.rg, L, rotation=(0, 90, 0) if along_x else (90, 0, 0))
    d = Vector(L / 2, 0, 0) if along_x else Vector(0, L / 2, 0)
    return Pos(x, y, zc) * (body + Pos(d) * Sphere(p.rg) + Pos(-d) * Sphere(p.rg))


def _ring_holes(n: int, r: float, d: float, z0: float, depth: float, phase: float = 0.0):
    return [
        Pos(r * math.cos(a), r * math.sin(a), z0) * Cylinder(d / 2, depth, align=MIN)
        for a in (phase + 2 * math.pi * k / n for k in range(n))
    ]


def make_parts(p: Params) -> dict[str, Part]:
    """Parts without balls at psi = alpha = 0. Names starting with '_' are bought."""
    c = p.slide_clr
    d2, D2, B2 = BEARINGS[p.bearing]
    zs = p.zs
    z4b, z5, z4a, z6, z7, z8 = (zs[k] for k in ("4b", "5", "4a", "6", "7", "8"))
    z_top = p.z_top
    ex = p.A  # eccentric offset at psi = 0 (+X)
    old = oldham_xy(p)
    R_s = p.R_s
    parts: dict[str, Part] = {}
    zb2 = z6 + p.t6 - p.bearing_lip_t - B2  # bearing 2 bottom (inserted from below)
    zb13 = p.zb13
    lip_r = D2 / 2 - p.bearing_shoulder

    # 1 input shaft: hex socket from below; eccentric cam in bearing 2; the flange on top of the
    # cam rests on the inner race (no axial drop) and stops under disc 8 (no lift-off).
    zs0 = p.z_bottom
    shaft = Pos(0, 0, zs0) * Cylinder(p.shaft_r, zb2 - zs0, align=MIN)
    shaft += Pos(ex, 0, zb2) * Cylinder(p.cam_r, B2, align=MIN)
    flange_r = d2 / 2 + (D2 - d2) / 8  # stays on the inner race
    shaft += Pos(ex, 0, zb2 + B2) * Cylinder(flange_r, z8 - c - zb2 - B2, align=MIN)
    hexs = extrude(RegularPolygon(p.hex_af / 2 + p.fit_clr, 6, major_radius=False), p.hex_depth)
    parts["1_input_shaft"] = shaft - Pos(0, 0, zs0) * hexs
    ring = Cylinder(D2 / 2, B2, align=MIN) - Cylinder(d2 / 2, B2, align=MIN)
    parts["_2_bearing_" + p.bearing] = Pos(ex, 0, zb2) * ring

    # 3 cover: bearing 13 in a boss below the plate carries the input shaft
    cover = _disc(0, p.R_cover, p.t3, lip_r)
    cover += _ring_holes(p.n_screws, R_s, 2 * p.ear_r, 0, p.t3)
    seat_r = D2 / 2 + p.fit_clr
    if zb13 < 0:
        cover += Pos(0, 0, zb13) * Cylinder(seat_r + p.bearing_boss_wall, -zb13, align=MIN)
    cover -= Pos(0, 0, zb13) * Cylinder(seat_r, B2, align=MIN)
    cover -= Pos(0, 0, zb13) * Cylinder(lip_r, p.t3 - zb13, align=MIN)
    parts["_13_bearing_" + p.bearing] = Pos(0, 0, zb13) * ring
    for x, y in old:  # Y slides on top face
        cover -= _slot(p, x, y, z4b, along_x=False)
    cover -= _ring_holes(p.n_screws, R_s, p.m3_clear_d, 0, p.t3)
    if p.cover_clamp < p.t3:  # counterbores so an available screw reaches the housing
        cover -= _ring_holes(p.n_screws, R_s, 2 * p.cbore_r, 0, p.t3 - p.cover_clamp)
    parts["3_end_cover"] = cover

    # 5 moves +-A in Y; its hole clears the shaft and the lower end of bearing 2
    cdisc = _disc(z5, p.R_out, p.t5, D2 / 2 + p.A + c)
    for x, y in old:
        cdisc -= _slot(p, x, y, z4b, along_x=False)
        cdisc -= _slot(p, x, y, z4a, along_x=True)
    parts["5_center_disc"] = cdisc

    hyp = _disc(z6, p.R_out, p.t6, lip_r)
    hyp -= Pos(0, 0, zb2) * Cylinder(seat_r, z6 + p.t6 - p.bearing_lip_t - zb2, align=MIN)
    for x, y in old:  # X slides on the left face
        hyp -= _slot(p, x, y, z4a, along_x=True)
    hyp -= _cycloid_groove(p, "hypo", z7, down=+1)
    parts["6_hypocycloid_disc"] = Pos(ex, 0, 0) * hyp

    # 8 output: hub through the lid bore, output face on top
    top = z_top + p.lid_t + p.hub_stub
    epi = _disc(z8, p.R_out, p.t8) - _cycloid_groove(p, "epi", z7, down=-1)
    hub_h = top - z8 - p.t8
    epi += Pos(0, 0, z8 + p.t8) * Cylinder(p.hub_r, hub_h, align=MIN)
    epi -= _ring_holes(p.n_out_holes, p.hub_r / 2, p.m3_tap_d, top - hub_h, hub_h)
    parts["8_epicycloid_disc_output"] = epi

    # 9 cage: round pockets; relative to the cage each ball circles with radius ~2A
    centers, pr = cage_pockets(p)
    r_in = min(abs(centers)) - pr - p.cage_ring
    r_out = max(abs(centers)) + pr + p.cage_ring
    cage = Pos(0, 0, z7) * (Cylinder(r_out, p.cage_t) - Cylinder(r_in, p.cage_t))
    cage -= [Pos(q.real, q.imag, z7) * Cylinder(pr, p.cage_t) for q in centers]
    parts["9_cage"] = cage

    t = p.t3
    housing = Pos(0, 0, t) * Cylinder(p.R_cover, z_top - t, align=MIN)
    housing += _ring_holes(p.n_screws, R_s, 2 * p.ear_r, t, z_top - t)  # screw ears
    housing -= Pos(0, 0, t) * Cylinder(p.R_h, z_top - t, align=MIN)
    hd_c, hd_l = p.hex_depth_10(p.cover_clamp), p.hex_depth_10(p.lid_t)
    hexagon = RegularPolygon(p.screw_d / 2, 6, major_radius=False)
    for a in range(p.n_screws):  # blind hex holes from both ends; screws cut into the corners
        ph = 2 * math.pi * a / p.n_screws
        for z0, hd in ((t, hd_c), (z_top - hd_l, hd_l)):
            housing -= Pos(R_s * math.cos(ph), R_s * math.sin(ph), z0) * extrude(hexagon, hd)
    parts["10_housing"] = housing
    lid = Pos(0, 0, z_top) * Cylinder(p.R_cover, p.lid_t, align=MIN)
    lid += _ring_holes(p.n_screws, R_s, 2 * p.ear_r, z_top, p.lid_t)
    bore = p.hub_r + p.lid_bore_clr
    lid += Pos(0, 0, z_top - p.lid_rib_h) * (  # thin rib: small contact ring for lift-off
        Cylinder(bore + p.lid_rib_w, p.lid_rib_h, align=MIN)
        - Cylinder(bore, p.lid_rib_h, align=MIN)
    )
    lid -= Pos(0, 0, z_top) * Cylinder(bore, p.lid_t, align=MIN)
    parts["11_lid"] = lid - _ring_holes(p.n_screws, R_s, p.m3_clear_d, z_top, p.lid_t)
    return parts


def cage_rotation(p: Params, psi: float) -> float:
    """Cage angle [rad]: mean angular travel of the balls 7."""
    return float(np.mean(np.angle(ball7_positions(p, psi) / ball7_positions(p, 0.0))))


def cage_pockets(p: Params, steps: int = 48) -> tuple[np.ndarray, float]:
    """Pocket centres (cage frame, psi=0) and radius covering each ball's path in the cage."""
    psis = np.linspace(0, 2 * np.pi, steps, endpoint=False)
    rel = np.array([ball7_positions(p, q) * np.exp(-1j * cage_rotation(p, q)) for q in psis])
    centers = rel.mean(axis=0)
    return centers, float(np.abs(rel - centers).max()) + p.rg + p.slide_clr


def oldham_xy(p: Params) -> list[tuple[float, float]]:
    ang = [math.pi / 4 + 2 * math.pi * k / p.n_oldham for k in range(p.n_oldham)]
    return [(p.Rb * math.cos(a), p.Rb * math.sin(a)) for a in ang]


def pose(p: Params, home: dict[str, Part], psi: float) -> dict[str, Part]:
    """Move the parts kinematically to input angle psi and add the balls."""
    e = p.A * np.exp(1j * psi)
    alpha = math.degrees(psi / p.ratio)
    b = ball7_positions(p, psi)
    cage_rot = math.degrees(cage_rotation(p, psi))
    move = {
        "1_input_shaft": Location((0, 0, 0), (0, 0, math.degrees(psi))),
        "_2_bearing_" + p.bearing: Pos(e.real - p.A, e.imag, 0),
        "5_center_disc": Pos(0, e.imag, 0),
        "6_hypocycloid_disc": Pos(e.real - p.A, e.imag, 0),
        "8_epicycloid_disc_output": Location((0, 0, 0), (0, 0, alpha)),
        "9_cage": Location((0, 0, 0), (0, 0, cage_rot)),
    }
    parts = {k: move[k] * v if k in move else v for k, v in home.items()}
    z7, z4b, z4a = p.zs["7"], p.zs["4b"], p.zs["4a"]
    parts["_7_reduction_balls"] = Compound(
        children=[Pos(c.real, c.imag, z7) * Sphere(p.rb) for c in b]
    )
    parts["_4_transmission_balls"] = Compound(
        children=[Pos(x + e.real / 2, y + e.imag, z4a) * Sphere(p.rb) for x, y in oldham_xy(p)]
        + [Pos(x, y + e.imag / 2, z4b) * Sphere(p.rb) for x, y in oldham_xy(p)]
    )
    return parts


def build_parts(p: Params) -> dict[str, Part]:
    return pose(p, make_parts(p), 0.0)


PRINT_FLIP = {"1_input_shaft", "3_end_cover", "11_lid"}  # turned over to avoid overhangs


def print_pose(name: str, part: Part) -> Part:
    """Rotate for support-free printing and drop the largest flat face onto z=0."""
    if name in PRINT_FLIP:
        part = part.rotate(Axis.X, 180)
    bb = part.bounding_box()
    return Pos(-bb.center().X, -bb.center().Y, -bb.min.Z) * part


def interference(parts: dict[str, Part], tol: float = 1e-3) -> list[tuple[str, str, float]]:
    """Pairs whose common volume exceeds tol (mm^3); ball/groove contact is surface-only."""
    names = list(parts)
    bad = []
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            bba, bbb = parts[a].bounding_box(), parts[b].bounding_box()
            if bba.min.Z >= bbb.max.Z or bbb.min.Z >= bba.max.Z:
                continue
            common = parts[a] & parts[b]
            v = common.volume if common else 0.0
            if v > tol:
                bad.append((a, b, v))
    return bad


def build(p: Params, out: Path) -> dict[str, Path]:
    """STEP for every part + assembly; STL (print pose) only for printed parts."""
    out.mkdir(parents=True, exist_ok=True)
    for old in [
        *out.glob("*.step"),
        *out.glob("*.stl"),
    ]:  # drop parts from earlier configs
        old.unlink()
    parts = build_parts(p)
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
    export_step(
        Compound(label="ball_reducer", children=list(parts.values())),
        str(out / "assembly.step"),
    )
    files["assembly.step"] = out / "assembly.step"
    asm_dir = out / "assembled"  # assembly-pose STLs for the viewer, bought parts included
    asm_dir.mkdir(exist_ok=True)
    for old in asm_dir.glob("*.stl"):
        old.unlink()
    for name, part in parts.items():
        f = asm_dir / f"{name.lstrip('_')}.stl"
        export_stl(part, str(f))
        files[f"assembled/{f.name}"] = f
    return files


if __name__ == "__main__":
    p = Params()
    d, D, B = BEARINGS[p.bearing]
    print(
        f"nb={p.nb} cage_t={p.cage_t:.2f} lid_t={p.lid_t:.2f} cover M3x{p.screw_len(p.cover_clamp)} lid M3x{p.screw_len(p.lid_t)} "
        f"t3/5/6/8={p.t3:.2f}/{p.t5:.2f}/{p.t6:.2f}/{p.t8:.2f} "
        f"height={p.height:.1f} "
        f"rho_min={min(min_curvature_radius(p, k) for k in ('epi', 'hypo')):.3f} "
        f"Rb={p.Rb:.3f} r0={p.r0:.4f} A={p.A:.4f} K={p.K} K2={p.K2} h={p.h:.3f} "
        f"groove_depth={p.groove_depth:.3f} R_out={p.R_out:.3f} outer_d={2 * p.R_cover:.2f} ear_d={2 * (p.R_s + p.ear_r):.2f} "
        f"shaft_r={p.shaft_r:.3f} bearing={p.bearing} {d}x{D}x{B} ratio={p.ratio}"
    )
    for k, v in wall_thicknesses(p).items():
        print(f"  wall {k}: {v:.2f}")
    for w in design_warnings(p) + [f"thin wall {k}: {v:.2f}" for k, v in thin_walls(p).items()]:
        print("WARNING:", w, file=sys.stderr)
    for f in build(p, Path("out")).values():
        print(f)
