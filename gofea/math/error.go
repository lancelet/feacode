package math

import "fmt"

const (
	errorIllegalIndex = iota + 1
	errorDimensionsIncompatible
)

type error_ int

func (e error_) Error() string {
	switch e {
	case errorIllegalIndex:
		return "Index out of bounds"
	case errorDimensionsIncompatible:
		return "Matrix dimensions are incompatible"
	}
	return fmt.Sprintf("Unknown error code %d", e)
}

var (
	ErrorIllegalIndex           error_ = error_(errorIllegalIndex)
	ErrorDimensionsIncompatible error_ = error_(errorDimensionsIncompatible)
)
