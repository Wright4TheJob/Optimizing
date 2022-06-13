import pytest
import FourBarOptimize as fb
import numpy as np
import math

def test_deltas_1():
  p1 = [0,0]
  p2 = [2, 1]
  assert fb.dx(p1, p2) == 2
  assert fb.dy(p1, p2) == 1

def test_deltas_2():
  p1 = [0,0]
  p2 = [-5, -9]
  assert fb.dy(p1, p2) == -9
  assert fb.dx(p1, p2) == -5

def test_deltas_trivial():
  p1 = [1,1]
  p2 = [1, 1]
  assert fb.dx(p1, p2) == 0
  assert fb.dy(p1, p2) == 0

def test_deltas_non_zero():
  p1 = [4,3]
  p2 = [5, 8]
  assert fb.dx(p1, p2) == 1
  assert fb.dy(p1, p2) == 5

def test_link_intersections_non_intersecting():
      base1 = [0,0]
      base2 = [10,0]
      with pytest.raises(ValueError) as e:
          fb.link_intersections(base1, base2, 1, 1)

def test_link_intersections_exact_1():
     base1 = [0,0]
     base2 = [10,0]
     assert fb.link_intersections(base1, base2, 5, 5) == [[5,0]]

def test_link_intersections_exact_non_equal():
     base1 = [0,0]
     base2 = [10,0]
     assert fb.link_intersections(base1, base2, 2, 8) == [[2,0]]

def test_link_intersections_exact_non_equal():
     base1 = [0,0]
     base2 = [2,0]
     assert fb.link_intersections(base1, base2, 2, 2)[0][0] == pytest.approx(1)
     assert fb.link_intersections(base1, base2, 2, 2)[0][1] == pytest.approx(math.sqrt(3))
     assert fb.link_intersections(base1, base2, 2, 2)[0][0] == pytest.approx(1)
     assert fb.link_intersections(base1, base2, 2, 2)[0][1] == pytest.approx(math.sqrt(3))

def test_law_of_cosines_equilateral():
    angle = np.pi/3
    assert fb.angle_from_triangle_lengths(1,1,1) == pytest.approx(angle)

def test_law_of_cosines_30_60_90():
    angle = 30*np.pi/180
    assert fb.angle_from_triangle_lengths(1,math.sqrt(3),2) == pytest.approx(angle)

def test_law_of_cosines_30_60_90_alt_order():
    angle = 30*np.pi/180
    assert fb.angle_from_triangle_lengths(1,2,math.sqrt(3)) == pytest.approx(angle)

def test_law_of_cosines_30_60_90_2():
    angle = 60*np.pi/180
    assert fb.angle_from_triangle_lengths(math.sqrt(3),1,2) == pytest.approx(angle)

def test_absolute_vector_angle_zero():
    assert fb.absolute_vector_angle([0,0], [1,0]) == 0

def test_absolute_vector_angle_q1():
    assert fb.absolute_vector_angle([0,0], [1,1]) == pytest.approx(np.pi/4)

def test_in_quadrent_1():
    assert fb.in_quadrent_1(1,1) == True
    assert fb.in_quadrent_1(0,0) == False
    assert fb.in_quadrent_1(1,-1) == False
    assert fb.in_quadrent_1(-1,1) == False
    assert fb.in_quadrent_1(-1,-1) == False

def test_in_quadrent_2():
    assert fb.in_quadrent_2(1,1) == False
    assert fb.in_quadrent_2(0,0) == False
    assert fb.in_quadrent_2(1,-1) == False
    assert fb.in_quadrent_2(-1,1) == True
    assert fb.in_quadrent_2(-1,-1) == False

def test_in_quadrent_3():
    assert fb.in_quadrent_3(1,1) == False
    assert fb.in_quadrent_3(0,0) == False
    assert fb.in_quadrent_3(1,-1) == False
    assert fb.in_quadrent_3(-1,1) == False
    assert fb.in_quadrent_3(-1,-1) == True

def test_in_quadrent_4():
    assert fb.in_quadrent_4(1,1) == False
    assert fb.in_quadrent_4(0,0) == False
    assert fb.in_quadrent_4(1,-1) == True
    assert fb.in_quadrent_4(-1,1) == False
    assert fb.in_quadrent_4(-1,-1) == False

def test_vector_angle_horizontal():
    assert fb.absolute_vector_angle([0,0], [1, 0]) == pytest.approx(0)

def test_vector_angle_horizontal_negative():
    assert fb.absolute_vector_angle([0,0], [-1, 0]) == pytest.approx(np.pi)

def test_vector_angle_vertical():
    assert fb.absolute_vector_angle([0,0], [0, 1]) == pytest.approx(np.pi/2)

def test_vector_angle_vertical_neg():
    assert fb.absolute_vector_angle([0,0], [0, -1]) == pytest.approx(3*np.pi/2)

def test_vector_angle_q1():
    assert fb.absolute_vector_angle([0,0], [1, 1]) == pytest.approx(np.pi/4)

def test_vector_angle_q2():
    assert fb.absolute_vector_angle([0,0], [-1, 1]) == pytest.approx(3*np.pi/4)

def test_vector_angle_q3():
    assert fb.absolute_vector_angle([0,0], [-1, -1]) == pytest.approx(5*np.pi/4)

def test_vector_angle_q4():
    assert fb.absolute_vector_angle([0,0], [1, -1]) == pytest.approx(7*np.pi/4)
