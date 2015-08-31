PRO INDX_COORD, q, q0

  COMMON PLUTO_GRID
  COMMON PLUTO_VAR
  PRINT,"> INDX_COORD, value = ",q0

  indx = WHERE(q eq q0)
  indx = ARRAY_INDICES(q, indx)
  PRINT,"i = ",string(indx[0],format='(i4)'),$
        "; j = ",string(indx[1],format='(i4)'),$
        "; k = ",string(indx[2],format='(i4)')
  PRINT,"x1,x2,x3 =",x1[indx[0]], x2[indx[1]], x3[indx[2]]


END

