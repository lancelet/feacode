package math

// Dense matrix.
type DenseMatrix struct {
	// Number of rows in the matrix.
	rows int
	// Number of columns in the matrix.
	cols int
	// Elements of the matrix, stored in row-major order.
	elements []float64
}

// Creates a dense matrix from a variable arguments list of its row-major
//  elements.
func FromRowMajorDenseV(rows, cols int,
	elements ...float64) (*DenseMatrix, error) {
	return FromRowMajorDense(rows, cols, elements)
}

// Creates a dense matrix from a slice of its row-major elements.  The
//  elements are copied.
func FromRowMajorDense(rows, cols int,
	elements []float64) (*DenseMatrix, error) {
	if rows*cols != len(elements) {
		return nil, ErrorDimensionsIncompatible
	}
	M := ZerosDense(rows, cols)
	copy(M.elements, elements)
	return M, nil
}

// Creates a dense matrix full of zeros.
func ZerosDense(rows, cols int) *DenseMatrix {
	return &DenseMatrix{rows, cols, make([]float64, rows*cols)}
}

// Sets an element in a dense matrix.
func (M *DenseMatrix) Set(r, c int, v float64) {
	z := r * M.cols
	M.elements[z : z+M.cols][c] = v
}

// Gets an element from a dense matrix.
func (M *DenseMatrix) Get(r, c int) float64 {
	z := r * M.cols
	return M.elements[z : z+M.cols][c]
}

// Checks if two dense matrices are equal up to an epsilon value.
func DenseApproxEq(A *DenseMatrix, B *DenseMatrix, eps float64) bool {
	if (A.rows != B.rows) || (A.cols != B.cols) {
		return false
	}
	for i := 0; i < (A.rows * A.cols); i++ {
		a := A.elements[i]
		b := B.elements[i]
		eq := (a+eps >= b) && (a-eps <= b)
		if !eq {
			return false
		}
	}
	return true
}

// Subtracts a dense matrix from another: C = A - B.
func SubDense(A *DenseMatrix, B *DenseMatrix) (*DenseMatrix, error) {
	if (A.rows != B.rows) || (A.cols != B.cols) {
		return nil, ErrorDimensionsIncompatible
	}
	C := ZerosDense(A.rows, A.cols)
	for i := 0; i < (A.rows * A.cols); i++ {
		C.elements[i] = A.elements[i] - B.elements[i]
	}
	return C, nil
}

// Multiplies two dense matrices: C = A * B.
func MulDense(A *DenseMatrix, B *DenseMatrix) (*DenseMatrix, error) {
	if A.cols != B.rows {
		return nil, ErrorDimensionsIncompatible
	}
	C := ZerosDense(A.rows, B.cols)
	for i := 0; i < A.rows; i++ {
		for j := 0; j < B.cols; j++ {
			sum := float64(0)
			for k := 0; k < A.cols; k++ {
				sum += A.elements[i*A.cols+k] * B.elements[k*B.cols+j]
			}
			C.elements[i*C.cols+j] = sum
		}
	}
	return C, nil
}
