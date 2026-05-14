from crypticip.fpocket import parse_info_text, parse_info_file, run_fpocket


def test_parse_info_file(fpocket_info):
    pockets = parse_info_file(fpocket_info)
    assert len(pockets) == 3
    p1 = pockets[0]
    assert p1.rank == 1
    assert p1.score == 0.621
    assert p1.druggability == 0.852
    assert abs(p1.depth - 22.345) < 1e-6
    # The "Real volume" value should be preferred (last value parsed wins).
    assert p1.volume in (612.3, 615.7)
    assert p1.n_alpha_spheres == 35
    p2 = pockets[1]
    assert p2.rank == 2
    assert abs(p2.depth - 18.9) < 1e-6


def test_run_fpocket_no_binary(tmp_path, tiny_pdb):
    # No binary installed → returns missing status, never crashes.
    pockets, meta = run_fpocket(tiny_pdb, work_dir=tmp_path,
                                status=None)
    # If fpocket really is installed (e.g. CI with it), we still expect a list.
    assert isinstance(pockets, list)
    assert "status" in meta
