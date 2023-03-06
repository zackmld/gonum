package mat

import (
	"gonum.org/v1/gonum/blas"
	"gonum.org/v1/gonum/blas/cblas128"
)

func (m *CDense) Add(a, b CMatrix) {
	ar, ac := a.Dims()
	br, bc := b.Dims()
	if ar != br || ac != bc {
		panic(ErrShape)
	}

	aU, aTrans, _ := untransposeExtractCmplx(a)
	bU, bTrans, _ := untransposeExtractCmplx(b)
	m.reuseAsNonZeroed(ar, ac)

	if arm, ok := a.(*CDense); ok {
		if brm, ok := b.(*CDense); ok {
			amat, bmat := arm.mat, brm.mat
			if m != aU {
				m.checkOverlap(amat)
			}
			if m != bU {
				m.checkOverlap(bmat)
			}
			for ja, jb, jm := 0, 0, 0; ja < ar*amat.Stride; ja, jb, jm = ja+amat.Stride, jb+bmat.Stride, jm+m.mat.Stride {
				for i, v := range amat.Data[ja : ja+ac] {
					m.mat.Data[i+jm] = v + bmat.Data[i+jb]
				}
			}
			return
		}
	}

	m.checkOverlapMatrix(aU)
	m.checkOverlapMatrix(bU)
	var restore func()
	if aTrans && m == aU {
		m, restore = m.isolatedWorkspace(aU)
		defer restore()
	} else if bTrans && m == bU {
		m, restore = m.isolatedWorkspace(bU)
		defer restore()
	}

	for r := 0; r < ar; r++ {
		for c := 0; c < ac; c++ {
			m.set(r, c, a.At(r, c)+b.At(r, c))
		}
	}
}

func (m *CDense) Sub(a, b CMatrix) {
	ar, ac := a.Dims()
	br, bc := b.Dims()
	if ar != br || ac != bc {
		panic(ErrShape)
	}

	aU, aTrans, _ := untransposeExtractCmplx(a)
	bU, bTrans, _ := untransposeExtractCmplx(b)
	m.reuseAsNonZeroed(ar, ac)

	if arm, ok := a.(*CDense); ok {
		if brm, ok := b.(*CDense); ok {
			amat, bmat := arm.mat, brm.mat
			if m != aU {
				m.checkOverlap(amat)
			}
			if m != bU {
				m.checkOverlap(bmat)
			}
			for ja, jb, jm := 0, 0, 0; ja < ar*amat.Stride; ja, jb, jm = ja+amat.Stride, jb+bmat.Stride, jm+m.mat.Stride {
				for i, v := range amat.Data[ja : ja+ac] {
					m.mat.Data[i+jm] = v - bmat.Data[i+jb]
				}
			}
			return
		}
	}

	m.checkOverlapMatrix(aU)
	m.checkOverlapMatrix(bU)
	var restore func()
	if aTrans && m == aU {
		m, restore = m.isolatedWorkspace(aU)
		defer restore()
	} else if bTrans && m == bU {
		m, restore = m.isolatedWorkspace(bU)
		defer restore()
	}

	for r := 0; r < ar; r++ {
		for c := 0; c < ac; c++ {
			m.set(r, c, a.At(r, c)-b.At(r, c))
		}
	}
}

func (m *CDense) Mul(a, b CMatrix) {

	ar, ac := a.Dims()
	br, bc := b.Dims()

	if ac != br {
		panic(ErrShape)
	}

	aU, aTrans, _ := untransposeExtractCmplx(a)
	bU, bTrans, _ := untransposeExtractCmplx(b)
	m.reuseAsNonZeroed(ar, bc)
	var restore func()
	if m == aU {
		m, restore = m.isolatedWorkspace(aU)
		defer restore()
	} else if m == bU {
		m, restore = m.isolatedWorkspace(bU)
		defer restore()
	}
	aT := blas.NoTrans
	if aTrans {
		aT = blas.Trans
	}
	bT := blas.NoTrans
	if bTrans {
		bT = blas.Trans
	}

	// Some of the cases do not have a transpose option, so create
	// temporary memory.
	// C = Aᵀ * B = (Bᵀ * A)ᵀ
	// Cᵀ = Bᵀ * A.
	if aU, ok := aU.(*CDense); ok {
		if restore == nil {
			m.checkOverlap(aU.mat)
		}
		switch bU := bU.(type) {
		case *CDense:
			if restore == nil {
				m.checkOverlap(bU.mat)
			}
			cblas128.Gemm(aT, bT, 1, aU.mat, bU.mat, 0, m.mat)
			return
		}
	}
}

func (m *CDense) Scale(f complex128, a CMatrix) {
	ar, ac := a.Dims()

	m.reuseAsNonZeroed(ar, ac)

	aU, aTrans, _ := untransposeExtractCmplx(a)
	if rm, ok := aU.(*CDense); ok {
		amat := rm.mat
		if m == aU || m.checkOverlap(amat) {
			var restore func()
			m, restore = m.isolatedWorkspace(a)
			defer restore()
		}
		if !aTrans {
			for ja, jm := 0, 0; ja < ar*amat.Stride; ja, jm = ja+amat.Stride, jm+m.mat.Stride {
				for i, v := range amat.Data[ja : ja+ac] {
					m.mat.Data[i+jm] = v * f
				}
			}
		} else {
			for ja, jm := 0, 0; ja < ac*amat.Stride; ja, jm = ja+amat.Stride, jm+1 {
				for i, v := range amat.Data[ja : ja+ar] {
					m.mat.Data[i*m.mat.Stride+jm] = v * f
				}
			}
		}
		return
	}

	m.checkOverlapMatrix(a)
	for r := 0; r < ar; r++ {
		for c := 0; c < ac; c++ {
			m.set(r, c, f*a.At(r, c))
		}
	}
}
