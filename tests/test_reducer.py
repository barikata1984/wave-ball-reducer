import math

import numpy as np
import pytest
from build123d import import_step, import_stl

from reducer import (
    Params,
    build,
    build_parts,
    crossing_counts,
    design_warnings,
    groove_xy,
    interference,
    make_parts,
    min_curvature_radius,
    pose,
    print_pose,
    search_K_K2,
    thin_walls,
)

PRINTED = {"1_input_shaft", "3_end_cover", "5_center_disc", "6_hypocycloid_disc",
           "8_epicycloid_disc_output", "9_cage", "10_housing", "11_lid"}  # fmt: skip


@pytest.fixture(scope="module")
def home():
    return make_parts(Params())


def test_derived_quantities():
    p = Params()
    assert p.K2 == pytest.approx(p.Rb * math.sin(math.pi / p.nb) / p.rb)
    assert p.A == pytest.approx(p.K * p.Rb / p.nb)
    assert p.h == pytest.approx(p.rb * (1 - math.sin(math.radians(p.beta_deg))))
    assert p.ratio == -10
    assert p.cage_t >= p.cage_min_t
    for clamp in (p.cover_clamp, p.lid_t):  # user owns M3x6, 8, 12, 16 only
        assert p.screw_len(clamp) in (6, 8, 12, 16) and p.screw_len(clamp) >= clamp + p.screw_engage


def test_ratio_is_half_of_epicycloid_waves():
    """Only alpha = -2 psi/(nb-1) keeps all ball crossings alive over a full input turn."""
    p = Params()
    assert min(crossing_counts(p, -(p.nb - 1) / 2)) >= 2 * p.nb
    for wrong in (-(p.nb - 1), p.nb - 1, (p.nb - 1) / 2, -(p.nb + 1) / 2):
        assert min(crossing_counts(p, wrong)) < 2 * p.nb, wrong


def test_defaults_come_from_search():
    p = Params()
    K, K2, d = search_K_K2(p)
    assert (K, K2) == (p.K, p.K2)
    assert d == pytest.approx(2 * p.R_cover) and d < 200
    assert p.A >= p.a_min
    for kind in ("epi", "hypo"):
        assert min_curvature_radius(p, kind) >= p.curv_margin * p.rg


def test_design_checks():
    assert not any("undercut" in w or "eccentricity" in w for w in design_warnings(Params()))
    w = design_warnings(Params(K=0.5, K2=1.3))
    assert any("undercut" in x for x in w) and any("K2" in x for x in w)
    assert any("eccentricity" in x for x in design_warnings(Params(K=0.1)))
    assert not any("Table I" in x for x in design_warnings(Params(nb=20, K=0.5)))


def test_min_wall():
    assert thin_walls(Params()) == {}
    assert "shaft around hex socket" in thin_walls(Params(hex_af=5.0))


def _assert_runs(p, home, psis):
    """Kinematic poses: printed parts never overlap, balls 7 sit in both grooves."""
    for psi in psis:
        parts = pose(p, home, psi)
        printed = {k: v for k, v in parts.items() if not k.startswith("_7")}
        assert interference(printed) == [], psi
        alpha, e = psi / p.ratio, p.A * np.exp(1j * psi)
        t = np.linspace(0, 2 * np.pi, 400000)
        c8 = np.exp(1j * alpha) * groove_xy(p, "epi", t)
        c6 = groove_xy(p, "hypo", t) + e
        for ball in parts["_7_reduction_balls"].solids():
            c = complex(ball.center().X, ball.center().Y)
            assert np.abs(c8 - c).min() < 2e-3 and np.abs(c6 - c).min() < 2e-3
        for disc in ("6_hypocycloid_disc", "8_epicycloid_disc_output"):
            for ball in parts["_7_reduction_balls"].solids()[:: max(1, p.nb // 7)]:
                assert ball.distance_to(parts[disc]) == pytest.approx(p.ball_clr, abs=0.03)
            assert parts["_7_reduction_balls"].distance_to(parts["9_cage"]) > 0
        gap = parts["11_lid"].distance_to(parts["8_epicycloid_disc_output"])
        assert gap == pytest.approx(min(p.lid_bore_clr, p.lid_axial_clr), abs=0.01)


def test_full_input_turn(home):
    _assert_runs(Params(), home, np.linspace(0, 2 * np.pi, 12, endpoint=False))


def test_print_pose_on_bed_and_fits_220(home):
    for name in PRINTED:
        part = print_pose(name, home[name])
        bb = part.bounding_box()
        assert bb.min.Z == pytest.approx(0, abs=1e-6)
        assert bb.size.X < 200 and bb.size.Y < 200 and bb.size.Z < 200
        assert any(abs(f.center().Z) < 1e-6 for f in part.faces())


def test_build_and_reload(tmp_path):
    files = build(Params(), tmp_path)
    stls = {f[:-4] for f in files if f.endswith(".stl") and "/" not in f}
    assert stls == PRINTED  # balls and bearings are bought: STEP only
    assert len([f for f in files if f.startswith("assembled/")]) == 12
    assert "assembly.step" in files
    for name, f in files.items():
        assert f.stat().st_size > 0
        if name.endswith(".step"):
            assert import_step(str(f)).volume > 0
        else:
            assert import_stl(str(f)).area > 0


def test_rebuild_other_nb(tmp_path):
    p = Params(nb=15, K2=2.6)
    files = build(p, tmp_path)
    assert len(import_step(str(files["7_reduction_balls.step"])).solids()) == 15
    assert interference(build_parts(p)) == []
    assert thin_walls(p) == {}
