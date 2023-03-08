// Copyright ©2018 The Gonum Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mat

import (
	"gonum.org/v1/gonum/blas/cblas128"
)

var (
	diagCDense *DiagCDense
	_          CMatrix      = diagCDense
	_          allMatrix    = diagCDense
	_          denseCMatrix = diagCDense
	_          CDiagonal    = diagCDense

	Cdiag CDiagonal
	_     CMatrix   = Cdiag
	_     CDiagonal = Cdiag
)

// Diagonal represents a diagonal matrix, that is a square matrix that only
// has non-zero terms on the diagonal.
type CDiagonal interface {
	CMatrix
	// Diag returns the number of rows/columns in the matrix.
	Diag() int

	// The following interfaces are included in the Diagonal
	// interface to allow the use of Diagonal types in
	// functions operating on these types.
}

// MutableDiagonal is a Diagonal matrix whose elements can be set.
type MutableCDiagonal interface {
	CDiagonal
	SetDiag(i int, v complex128)
}

// DiagCDense represents a diagonal matrix in dense storage format.
type DiagCDense struct {
	mat cblas128.Vector
}

// NewDiagDense creates a new Diagonal matrix with n rows and n columns.
// The length of data must be n or data must be nil, otherwise NewDiagDense
// will panic. NewDiagDense will panic if n is zero.
func NewDiagCDense(n int, data []complex128) *DiagCDense {
	if n <= 0 {
		if n == 0 {
			panic(ErrZeroLength)
		}
		panic("mat: negative dimension")
	}
	if data == nil {
		data = make([]complex128, n)
	}
	if len(data) != n {
		panic(ErrShape)
	}
	return &DiagCDense{
		mat: cblas128.Vector{N: n, Data: data, Inc: 1},
	}
}

// Diag returns the dimension of the receiver.
func (d *DiagCDense) Diag() int {
	return d.mat.N
}

// Dims returns the dimensions of the matrix.
func (d *DiagCDense) Dims() (r, c int) {
	return d.mat.N, d.mat.N
}

// T returns the transpose of the matrix.
func (d *DiagCDense) T() CMatrix {
	return d
}

func (d *DiagCDense) H() CMatrix {
	return ConjTranspose{d}
}

// Bandwidth returns the upper and lower bandwidths of the matrix.
// These values are always zero for diagonal matrices.
func (d *DiagCDense) Bandwidth() (kl, ku int) {
	return 0, 0
}

// Reset empties the matrix so that it can be reused as the
// receiver of a dimensionally restricted operation.
//
// Reset should not be used when the matrix shares backing data.
// See the Reseter interface for more information.
func (d *DiagCDense) Reset() {
	// No change of Inc or n to 0 may be
	// made unless both are set to 0.
	d.mat.Inc = 0
	d.mat.N = 0
	d.mat.Data = d.mat.Data[:0]
}

// Zero sets all of the matrix elements to zero.
func (d *DiagCDense) Zero() {
	for i := 0; i < d.mat.N; i++ {
		d.mat.Data[d.mat.Inc*i] = 0
	}
}

// DiagView returns the diagonal as a matrix backed by the original data.
func (d *DiagCDense) DiagView() CDiagonal {
	return d
}

// DiagFrom copies the diagonal of m into the receiver. The receiver must
// be min(r, c) long or empty, otherwise DiagFrom will panic.
func (d *DiagCDense) DiagFrom(m CMatrix) {
	n := min(m.Dims())
	d.reuseAsNonZeroed(n)

	var vec cblas128.Vector
	switch r := m.(type) {
	case *DiagCDense:
		vec = r.mat
	default:
		for i := 0; i < n; i++ {
			d.setDiag(i, m.At(i, i))
		}
		return
	}
	cblas128.Copy(vec, d.mat)
}

// RawBand returns the underlying data used by the receiver represented
// as a cblas128.Band.
// Changes to elements in the receiver following the call will be reflected
// in returned cblas128.Band.
func (d *DiagCDense) RawBand() cblas128.Band {
	return cblas128.Band{
		Rows:   d.mat.N,
		Cols:   d.mat.N,
		KL:     0,
		KU:     0,
		Stride: d.mat.Inc,
		Data:   d.mat.Data,
	}
}

// reuseAsNonZeroed resizes an empty diagonal to a r×r diagonal,
// or checks that a non-empty matrix is r×r.
func (d *DiagCDense) reuseAsNonZeroed(r int) {
	if r == 0 {
		panic(ErrZeroLength)
	}
	if d.IsEmpty() {
		d.mat = cblas128.Vector{
			Inc:  1,
			Data: useC(d.mat.Data, r),
		}
		d.mat.N = r
		return
	}
	if r != d.mat.N {
		panic(ErrShape)
	}
}

// IsEmpty returns whether the receiver is empty. Empty matrices can be the
// receiver for size-restricted operations. The receiver can be emptied using
// Reset.
func (d *DiagCDense) IsEmpty() bool {
	// It must be the case that d.Dims() returns
	// zeros in this case. See comment in Reset().
	return d.mat.Inc == 0
}

// Trace returns the trace of the matrix.
//
// Trace will panic with ErrZeroLength if the matrix has zero size.
/*func (d *DiagCDense) Trace() complex128 {
	if d.IsEmpty() {
		panic(ErrZeroLength)
	}
	rb := d.RawBand()
	var tr complex128
	for i := 0; i < rb.Rows; i++ {
		tr += rb.Data[rb.KL+i*rb.Stride]
	}
	return tr
}*/

// Norm returns the specified norm of the receiver. Valid norms are:
//
//	1 or Inf - The maximum diagonal element magnitude
//	2 - The Frobenius norm, the square root of the sum of the squares of
//	    the diagonal elements
//
// Norm will panic with ErrNormOrder if an illegal norm is specified and with
// ErrZeroLength if the receiver has zero size.
func (d *DiagCDense) Norm(norm complex128) complex128 {
	/*if d.IsEmpty() {
		panic(ErrZeroLength)
	}
	switch norm {
	default:
		panic(ErrNormOrder)
	case 1, math.Inf(1):
		imax := cblas128.Iamax(d.mat)
		return cmplx.Abs(d.at(imax, imax))
	case 2:
		return cblas128.Nrm2(d.mat)
	}*/
	return 0
}
