;================================ TIDY ===============================
;This routine sets everything back to idl default values to clean up
;at the end of a program.

PRO tidy

set_plot, 'X'
!P.FONT=-1
!P.multi=0
!P.THICK=1.
!X.THICK=1.
!Y.THICK=1.
!Z.THICK=1.
!P.CHARTHICK=1.
!P.CHARSIZE=1.
!X.STYLE=0
!Y.STYLE=0

END
