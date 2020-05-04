import pytest
from pytest import approx
from hephaestus import thermocouples

examples = [
    ("B", 0.291279, 250.0),
    ("B", 1.974546, 630.0),
    ("B", 1.980771, 631.0),
    ("B", 5.779517, 1100.0),
    ("B", 13.820279, 1820.0),
    ("B", 0.006197, 60.0),
    ("C", 0.000000, 0.0),
    ("C", 1.451164, 100.0),
    ("C", 11.778834, 660.0),
    ("C", 20.066477, 1100.0),
    ("C", 29.402817, 1680.0),
    ("C", 37.061295, 2315.0),
    ("E", -9.834951, -270.0),
    ("E", -8.824581, -200.0),
    ("E", -5.237184, -100.0),
    ("E", 0.000000, 0.0),
    ("E", 6.318930, 100.0),
    ("E", 49.917224, 660.0),
    ("E", 76.372826, 1000.0),
    ("J", -8.095380, -210.0),
    ("J", -4.632524, -100.0),
    ("J", 0.000000, 0.0),
    ("J", 36.675405, 660.0),
    ("J", 43.559497, 770.0),
    ("J", 69.553180, 1200.0),
    ("K", -6.457738, -270.0),
    ("K", -5.891404, -200.0),
    ("K", -3.553631, -100.0),
    ("K", 0.000000, 0.0),
    ("K", 5.206093, 127.0),
    ("K", 27.447068, 660.0),
    ("K", 54.886364, 1372.0),
    ("N", -4.3451354, -270.0),
    ("N", -3.990376, -200.0),
    ("N", -2.406811, -100.0),
    ("N", 0.000000, 0.0),
    ("N", 2.774124, 100.0),
    ("N", 22.957828, 660.0),
    ("N", 47.512772, 1300.0),
    ("R", -0.226465, -50.0),
    ("R", -0.123155, -25.0),
    ("R", 0.000000, 0.0),
    ("R", 6.273327, 660.0),
    ("R", 11.849642, 1100.0),
    ("R", 21.101477, 1768.0),
    ("S", -0.235555, -50.0),
    ("S", -0.126831, -25.0),
    ("S", 0.000000, 0.0),
    ("S", 5.856769, 660.0),
    ("S", 10.756545, 1100.0),
    ("S", 18.692510, 1768.0),
    ("T", -6.257505, -270.0),
    ("T", -5.602961, -200.0),
    ("T", -3.378582, -100.0),
    ("T", 0.000000, 0.0),
    ("T", 4.278519, 100.0),
    ("T", 9.288102, 200.0),
    ("T", 20.871970, 400.0),
]


@pytest.mark.parametrize("ttype, emf, t90", examples)
def test_emfr(ttype, emf, t90):
    tc = thermocouples[ttype]
    assert tc.emfr(t90) == approx(emf, abs=1e-6)


@pytest.mark.parametrize("ttype, emf, t90", examples)
def test_t90r(ttype, emf, t90):
    tc = thermocouples[ttype]
    assert tc.t90r(emf) == approx(t90, abs=1e-2)
