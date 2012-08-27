package math

import (
	. "launchpad.net/gocheck"
	"testing"
)

// Hook up gocheck into the gotest runner.
func Test(t *testing.T) { TestingT(t) }

type MatrixTestSuite struct{}

var _ = Suite(&MatrixTestSuite{})

func (s *MatrixTestSuite) TestFromRowMajorDenseV(c *C) {
	expected := &DenseMatrix{2, 2, []float64{1, 2, 3, 4}}
	m, merr := FromRowMajorDenseV(2, 2, 1, 2, 3, 4)
	c.Check(merr, IsNil)
	c.Check(DenseApproxEq(expected, m, 1.0E-12), Equals, true)
	q, qerr := FromRowMajorDenseV(2, 2, 1, 2, 3, 4, 5)
	c.Check(q, IsNil)
	c.Check(qerr, Equals, ErrorDimensionsIncompatible)
}

func (s *MatrixTestSuite) TestFromRowMajorDense(c *C) {
	els := []float64{1, 2, 3, 4, 5, 6}
	m, err := FromRowMajorDense(2, 3, els)
	c.Check(err, IsNil)
	c.Check(m.rows, Equals, 2)
	c.Check(m.cols, Equals, 3)
	c.Check(m.Get(0, 0), Equals, float64(1))
	c.Check(m.Get(0, 1), Equals, float64(2))
	c.Check(m.Get(0, 2), Equals, float64(3))
	c.Check(m.Get(1, 0), Equals, float64(4))
	c.Check(m.Get(1, 1), Equals, float64(5))
	c.Check(m.Get(1, 2), Equals, float64(6))
	q, qerr := FromRowMajorDense(2, 2, els)
	c.Check(q, IsNil)
	c.Check(qerr, Equals, ErrorDimensionsIncompatible)
}

func (s *MatrixTestSuite) TestZerosDense(c *C) {
	const nRows = 4
	const nCols = 3
	z := ZerosDense(nRows, nCols)
	c.Check(z.rows, Equals, nRows)
	c.Check(z.cols, Equals, nCols)
	for row := 0; row < nRows; row++ {
		for col := 0; col < nCols; col++ {
			c.Check(z.Get(row, col), Equals, float64(0))
		}
	}
}

func (s *MatrixTestSuite) TestSet(c *C) {
	const nRows = 7
	const nCols = 8
	const setRow = 2
	const setCol = 1
	const setValue = 42.0
	m := ZerosDense(nRows, nCols)
	m.Set(setRow, setCol, setValue)
	for row := 0; row < nRows; row++ {
		for col := 0; col < nCols; col++ {
			expected := float64(0)
			if (row == setRow) && (col == setCol) {
				expected = setValue
			}
			c.Check(m.Get(row, col), Equals, expected)
		}
	}
}

func (s *MatrixTestSuite) TestDenseApproxEq(c *C) {
	a, _ := FromRowMajorDenseV(2, 2, 1, 2, 3, 4)
	b, _ := FromRowMajorDenseV(2, 2, 1, 2.2, 3, 4)
	c.Check(DenseApproxEq(a, b, 0.2001), Equals, true)
	c.Check(DenseApproxEq(a, b, 0.1999), Equals, false)
	q := ZerosDense(3, 2)
	c.Check(DenseApproxEq(a, q, 1000), Equals, false)
}

func (s *MatrixTestSuite) TestSubDense(c *C) {
	a, erra := FromRowMajorDenseV(2, 2, 1, 2, 3, 4)
	b, errb := FromRowMajorDenseV(2, 2, 7, 2, 1, 9)
	expected, errExpected := FromRowMajorDenseV(2, 2, -6, 0, 2, -5)
	c.Check(erra, IsNil)
	c.Check(errb, IsNil)
	c.Check(errExpected, IsNil)
	m, err := SubDense(a, b)
	c.Check(err, IsNil)
	c.Check(DenseApproxEq(m, expected, 1.0E-12), Equals, true)
	d := ZerosDense(3, 2)
	q, errd := SubDense(a, d)
	c.Check(q, IsNil)
	c.Check(errd, Equals, ErrorDimensionsIncompatible)
}

func (s *MatrixTestSuite) TestMulDense(c *C) {
	a, erra := FromRowMajorDenseV(2, 2, 1, 2, 3, 4)
	b, errb := FromRowMajorDenseV(2, 2, 7, 2, 1, 9)
	expected, errExpected := FromRowMajorDenseV(2, 2, 9, 20, 25, 42)
	c.Check(erra, Equals, nil)
	c.Check(errb, Equals, nil)
	c.Check(errExpected, Equals, nil)
	m, err := MulDense(a, b)
	c.Check(err, Equals, nil)
	c.Check(DenseApproxEq(m, expected, 1.0E-12), Equals, true)
	d := ZerosDense(3, 2)
	q, errd := MulDense(a, d)
	c.Check(q, IsNil)
	c.Check(errd, Equals, ErrorDimensionsIncompatible)
}
