"""Parametric part builders shared by the reducers (build123d).

Every function takes its dimensions as arguments and returns build123d solids in the assembled
pose: axis = Z, millimetres, angles in degrees unless the name says otherwise. Functions whose
name ends in ``_cut`` return a cutter to subtract; the others return material. Nothing here
knows about a particular reducer, so stepped.py, planetary.py or reducer.py can all call it.

FDM conventions assumed throughout: overhangs no flatter than 45 deg, flat ceilings only as
short bridges, clearances radial unless stated.
"""

from __future__ import annotations

import math

from bd_warehouse.thread import Thread
from build123d import (
    Align,
    Axis,
    Box,
    Cone,
    Cylinder,
    Face,
    Part,
    Polyline,
    Pos,
    RegularPolygon,
    Rot,
    Solid,
    Wire,
    extrude,
    revolve,
)

MIN = (None, None, Align.MIN)  # as reducer.MIN: x, y as built, z from 0
# The bd_warehouse ridge is a B-spline sweep; where its root corner lies exactly on the base
# cylinder, the fuse leaves slivers ~0.005 wide. Bury the corner: the base reaches this far
# past the thread root.
ROOT_BURY = 0.02


# ---------------------------------------------------------------- small primitives
def sector(r0: float, r1: float, deg0: float, deg1: float, z0: float, h: float) -> Part:
    """Annular sector r0..r1, angles deg0..deg1, z0..z0 + h (true arcs, so cuts line up)."""
    rect = Face(Wire(Polyline((r0, 0, z0), (r1, 0, z0), (r1, 0, z0 + h), (r0, 0, z0 + h),
                              close=True).edges()))  # fmt: skip
    return revolve(rect, Axis.Z, deg1 - deg0).rotate(Axis.Z, deg0)


def hex_prism(af: float, z0: float, h: float, x: float = 0.0, y: float = 0.0) -> Part:
    """Hexagonal prism, across-flats af, from z0 up by h, centred on (x, y)."""
    return Pos(x, y, z0) * extrude(RegularPolygon(af / 2, 6, major_radius=False), h)


def ring_holes(n: int, r: float, d: float, z0: float, depth: float, phase_deg: float = 0.0):
    """n cylinders d x depth from z0 on a circle of radius r (cutters)."""
    return [Pos(r * math.cos(a), r * math.sin(a), z0) * Cylinder(d / 2, depth, align=MIN)
            for a in (math.radians(phase_deg) + 2 * math.pi * k / n for k in range(n))]  # fmt: skip


def bearing_dummy(d: float, D: float, B: float, z0: float) -> Part:
    """Bought ball bearing as a plain ring d..D, z0..z0 + B (races and balls not modelled)."""
    return Pos(0, 0, z0) * (Cylinder(D / 2, B, align=MIN) - Cylinder(d / 2, B, align=MIN))


# ---------------------------------------------------------------- threads and screw-on cap
def thread_45(
    r_root: float, depth: float, pitch: float, z0: float, length: float, internal: bool = False
) -> Solid:
    """Right-hand trapezoid thread ridge with 45 deg flanks (90 deg included), via bd_warehouse.

    The ridge is 2 depth + 0.2 wide at its root and 0.2 at its crest (flats 0.2, so it prints
    without support). External: root on radius r_root, crest at r_root + depth. Internal: root
    on r_root + depth, crest at r_root (r_root is then the crest/minor radius). Both ends are cut
    square at z0 and z0 + length. At angle 0 the ridge centres lie at z0 + pitch/2 + k pitch, and
    the helix rises with +angle (right hand). The root reaches 0.2 into the base for the fuse.
    """
    crest, root = (r_root + depth, r_root) if not internal else (r_root, r_root + depth)
    return Pos(0, 0, z0) * Thread(apex_radius=crest, apex_width=0.2, root_radius=root,
                                  root_width=2 * depth + 0.2, pitch=pitch, length=length,
                                  end_finishes=("square", "square"))  # fmt: skip


def housing_thread(r_root: float, depth: float, pitch: float, z_top: float, turns: float) -> Solid:
    """External ridge for the housing top band, ending flush at z_top (the cap's seat face).

    Fuse onto a cylinder of radius r_root + ROOT_BURY. Length (turns) x pitch below z_top; make it at least
    one turn longer than the cap thread so the cap ridge never meets the start of this ridge.
    Ridge centres at angle 0 lie at z_top + (k + 1/2) pitch, i.e. midway between the cap's.
    """
    return thread_45(r_root, depth, pitch, z_top - turns * pitch, turns * pitch)


def screw_cap(
    r_thread_root: float,
    thread_depth: float,
    pitch: float,
    thread_clr: float,
    wall: float,
    z_seat: float,
    plate_t: float,
    turns: float,
    rib_r: float = 0.0,
    rib_h: float = 0.0,
) -> Part:
    """Cap screwed onto an external thread made by housing_thread(r_thread_root, ...).

    Returns the cap seated: the plate underside lies on z_seat (the housing top face), which
    fixes its height whatever the thread does. The skirt bore is the housing crest radius +
    thread_clr; the internal ridge (crest at r_thread_root + thread_clr) sits midway in the
    housing groove, so with 45 deg flanks there is thread_clr of axial play each way (flank
    normal gap thread_clr sin 45). Outer radius = crest + thread_clr + wall. The optional rib
    (disc rib_r x rib_h under the plate) reaches into the housing bore, e.g. down to a carrier.
    Assumes it prints plate-down (upside down): skirt, ridges and rib then need no support.
    """
    r_in = r_thread_root + thread_depth + thread_clr
    r_out = r_in + wall
    skirt_h = (turns + 0.5) * pitch
    cap = Pos(0, 0, z_seat) * Cylinder(r_out, plate_t, align=MIN)
    cap += Pos(0, 0, z_seat - skirt_h) * (Cylinder(r_out, skirt_h, align=MIN)
                                          - Cylinder(r_in - ROOT_BURY, skirt_h, align=MIN))  # fmt: skip
    cap += thread_45(r_thread_root + thread_clr, thread_depth, pitch, z_seat - skirt_h, skirt_h,
                     internal=True)  # fmt: skip
    if rib_h > 0:
        cap += Pos(0, 0, z_seat - rib_h) * Cylinder(rib_r, rib_h, align=MIN)
    return cap


def cap_outer_r(r_thread_root: float, thread_depth: float, thread_clr: float, wall: float) -> float:
    """Outer radius of screw_cap for these arguments."""
    return r_thread_root + thread_depth + thread_clr + wall


def ratchet_arm(
    cap: Part,
    r_out: float,
    width: float,
    slot: float,
    end_deg: float,
    span_deg: float,
    z_face: float,
    plate_t: float,
    below_h: float,
    tooth_deg: float,
    tooth_r: tuple[float, float],
    tooth_h: float,
    tooth_flat: float,
) -> Part:
    """Cut a springy arc arm into the rim of a plate and hang a ratchet tooth under its free end.

    The arm is the plate ring r_out - width .. r_out between end_deg (free end) and
    end_deg + span_deg (root), cut free by an arc slot inside it and a radial slot at its end;
    everything of `cap` below the arm down to z_face - below_h (a skirt) is removed, so the arm
    can be lifted at the end to release the lock. The tooth (tooth_r radial extent, tooth_h
    deep) hangs under z_face; its vertical face is at tooth_deg and faces +theta (blocks the
    loosening turn of a right-hand thread), followed by a flat and a 45 deg ramp towards -theta.
    Returns the modified cap.
    """
    a0, a1 = end_deg, end_deg + span_deg
    r_in = r_out - width
    rs = r_in - slot
    ds = math.degrees(slot / r_out)
    cap -= sector(rs, r_in, a0 - ds, a1, z_face, plate_t)
    cap -= sector(rs, r_out + 1, a0 - ds, a0, z_face, plate_t)
    cap -= sector(r_in, r_out + 1, a0 - ds, a1, z_face - below_h, below_h)
    th, tf = tooth_h, tooth_flat
    prof = Face(Wire(Polyline((0, 0, 0), (0, 0, -th), (0, -tf, -th), (0, -tf - th, 0),
                              close=True).edges()))  # fmt: skip
    r0, r1 = tooth_r
    tooth = Pos(r0, 0, z_face) * extrude(prof, r1 - r0, (1, 0, 0))
    return cap + tooth.rotate(Axis.Z, tooth_deg)


def ratchet_notches_cut(
    r0: float,
    r1: float,
    z_top: float,
    n: int,
    tooth_deg: float,
    r_mid: float,
    tooth_h: float,
    tooth_flat: float,
    play: float,
) -> list:
    """n notches in a face at z_top (cutters), one in front of the tooth of ratchet_arm.

    Each notch spans the tooth (flat + ramp) plus `play` on both sides at radius r_mid and is
    tooth_h + 0.2 deep; its +theta wall is what the tooth's vertical face hits when loosening.
    """
    w0 = math.degrees(play / r_mid)
    w1 = math.degrees((tooth_flat + tooth_h + play) / r_mid)
    dn = tooth_h + 0.2
    return [sector(r0, r1, c - w1, c + w0, z_top - dn, dn + 1)
            for c in (tooth_deg + 360 * k / n for k in range(n))]  # fmt: skip


def grip_grooves_cut(r: float, z0: float, h: float, n: int, w: float, d: float) -> list:
    """n vertical finger grooves w wide, d deep into a rim of radius r, z0..z0 + h (cutters)."""
    return [Rot(0, 0, 360 * k / n) * Pos(r - d, -w / 2, z0) * Box(d + 1, w, h, align=MIN)
            for k in range(n)]  # fmt: skip


# ---------------------------------------------------------------- split housing
def split_housing(
    housing: Part,
    z_split: float,
    r_out: float,
    spigot_r0: float,
    spigot_w: float,
    spigot_h: float,
    spigot_in_clr: float,
    slide_clr: float,
    n_bolts: int,
    bolt_r: float,
    bolt_phase_deg: float,
    clamp: float,
    cbore_top: float,
    screw_d: float,
    clear_d: float,
    head_d: float,
    screw_len: float,
    hex_extra: float = 1.0,
) -> tuple[Part, Part, Part]:
    """Split a one-piece housing at z_split into (upper, lower, spigot ring) joined by axial
    cap screws.

    Both split faces stay flat (each is a bed face: the upper prints upright, the lower upside
    down), so the spigot is a loose ring spigot_r0 .. spigot_r0 + spigot_w, 2 spigot_h high,
    half in a groove in each part. Each groove's outer wall has slide_clr (these walls centre the
    parts), its inner wall spigot_in_clr (no contact), its floor 0.2 past the ring and its roof a
    45 deg ridge, so both parts print on the split face without support. The ring must not
    cross the bolt holes (put it outside or inside the bolt circle).
    Bolts: n_bolts screws on radius bolt_r, put in from above; clearance hole clear_d through
    the upper's slab z_split .. z_split + clamp, counterbore head_d + 2 slide_clr from there up
    to cbore_top (must be open above for the key), hex hole (across flats screw_d, the screw
    cuts its own thread) screw_len - clamp + hex_extra deep into the lower.
    r_out only sizes the cutting boxes (>= the housing radius).
    """
    big = 2 * r_out + 2
    bb = housing.bounding_box()
    upper = housing & (Pos(-big / 2, -big / 2, z_split) * Box(big, big, bb.max.Z - z_split,
                                                              align=MIN))  # fmt: skip
    lower = housing & (Pos(-big / 2, -big / 2, bb.min.Z) * Box(big, big, z_split - bb.min.Z,
                                                               align=MIN))  # fmt: skip
    r0, r1 = spigot_r0, spigot_r0 + spigot_w
    ring = Pos(0, 0, z_split - spigot_h) * (Cylinder(r1, 2 * spigot_h, align=MIN)
                                            - Cylinder(r0, 2 * spigot_h, align=MIN))  # fmt: skip
    g0, g1 = r0 - spigot_in_clr, r1 + slide_clr
    zw = spigot_h + 0.2
    zr = zw + (g1 - g0) / 2

    def groove(sz: int) -> Part:  # sz = +1 into the upper, -1 into the lower
        prof = Polyline((g0, 0, z_split), (g1, 0, z_split), (g1, 0, z_split + sz * zw),
                        ((g0 + g1) / 2, 0, z_split + sz * zr), (g0, 0, z_split + sz * zw),
                        close=True)  # fmt: skip
        return revolve(Face(Wire(prof.edges())), Axis.Z, 360)

    upper -= groove(1)
    lower -= groove(-1)
    upper -= ring_holes(n_bolts, bolt_r, clear_d, z_split, clamp, bolt_phase_deg)
    upper -= ring_holes(n_bolts, bolt_r, head_d + 2 * slide_clr, z_split + clamp,
                        cbore_top - z_split - clamp, bolt_phase_deg)  # fmt: skip
    hd = screw_len - clamp + hex_extra
    for k in range(n_bolts):
        a = math.radians(bolt_phase_deg) + 2 * math.pi * k / n_bolts
        lower -= hex_prism(screw_d, z_split - hd, hd, bolt_r * math.cos(a), bolt_r * math.sin(a))
    return upper, lower, ring


# ---------------------------------------------------------------- sun drive
def sun_adapter(
    r: float,
    t: float,
    seat_r: float,
    seat_top: float,
    hex_af: float,
    hex_top: float,
    n_screws: int,
    pcd: float,
    clear_d: float,
    phase_deg: float = 0.0,
) -> Part:
    """Adapter screwed to a rotor face at z = 0: disc r x t, round seat, hex drive.

    The round seat (radius seat_r, top at seat_top) carries the sun's bottom face, so the hex
    (across flats hex_af, up to hex_top) only turns the sun. Screw holes clear_d on pcd from
    the disc top (heads sit on the disc; the seat must stay inside the head circle).
    Prints disc-down without support.
    """
    adapter = Cylinder(r, t, align=MIN)
    adapter += Pos(0, 0, t) * Cylinder(seat_r, seat_top - t, align=MIN)
    adapter += hex_prism(hex_af, t, hex_top - t)
    return adapter - ring_holes(n_screws, pcd / 2, clear_d, 0, t, phase_deg)


def hex_socket_cut(hex_af: float, hex_clr: float, z0: float, depth: float, chamfer: float) -> Part:
    """Hex socket (across flats hex_af + 2 hex_clr) from z0 up by depth, 45 deg mouth chamfer.

    The chamfer (chamfer wide at z0) takes the first-layer squish of a socket printed mouth
    down. The roof is a hex bridge of hex_af + 2 hex_clr across flats.
    """
    sock = hex_af + 2 * hex_clr
    mouth = RegularPolygon(sock / 2 + chamfer, 6, major_radius=False)
    return hex_prism(sock, z0, depth) + Pos(0, 0, z0) * extrude(mouth, chamfer, taper=45)


# ---------------------------------------------------------------- bearings and pins
def bearing_seat_cut(D: float, B: float, fit_clr: float, z_mouth: float, z_floor: float) -> Part:
    """Pocket for a bearing (outer diameter D, width B) opening upward at z_mouth.

    Seat radius D/2 + fit_clr over B, then a 45 deg cone narrowing down to z_floor, on which the
    outer race's edge rests (z_floor < z_mouth - B). Printed mouth down, the cone needs no
    support; the flat ceiling left at z_floor (outer radius D/2 + fit_clr - (z_mouth - B -
    z_floor)) is a bridge around whatever bore the caller cuts through it.
    """
    r = D / 2 + fit_clr
    h = z_mouth - B - z_floor
    return (Pos(0, 0, z_mouth - B) * Cylinder(r, B, align=MIN)
            + Pos(0, 0, z_floor) * Cone(r - h, r, h, align=MIN))  # fmt: skip


def printed_pin(
    r: float,
    z_bottom: float,
    z_top: float,
    chamfer: float = 0.5,
    collar_r: float = 0.0,
    collar_h: float = 0.0,
) -> Part:
    """Printed pin r from z_bottom to z_top (fuse at z_top onto a plate), chamfered free end.

    collar_r x collar_h at the top end is a shoulder, e.g. to hold a bearing inner race (keep
    collar_r inside the inner-race land). Printed with the free end up, it needs no support.
    """
    pin = Pos(0, 0, z_bottom + chamfer) * Cylinder(r, z_top - z_bottom - chamfer, align=MIN)
    pin += Pos(0, 0, z_bottom) * Cone(r - chamfer, r, chamfer, align=MIN)
    if collar_h > 0:
        pin += Pos(0, 0, z_top - collar_h) * Cylinder(collar_r, collar_h, align=MIN)
    return pin
