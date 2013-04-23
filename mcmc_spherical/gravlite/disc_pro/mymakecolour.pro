;============================ MYMAKECOLOUR ============================
;This routine geneartes color table such that:
;color=0 is black
;color=1 is white
;color=2 is red
;color=3 is green
;color=4 is blue
;color=5 is yellow
;color=6 is magenta
;color=7 is cyan
;color=8 is pink
;color=9 is violet
;color=10 is lime
;color=11 is brown
;color=12 is sea blue
;color=13 is red pink
;color=14 is purple

PRO mymakecolour

;Seven colours for now, but can easily add more!
ncolours=17

;Set up colour table. Use RGB values:
RGB=fltarr(ncolours+1,3)
RGB(0,*)=[0,0,0]                ;Black
RGB(1,*)=[255,255,255]          ;White
RGB(2,*)=[255,0,0]              ;Red
RGB(3,*)=[0,255,0]              ;Green
RGB(4,*)=[0,0,255]              ;Blue
RGB(5,*)=[255,255,0]            ;Yellow
RGB(6,*)=[255,0,255]            ;Magenta
RGB(7,*)=[0,255,255]            ;Cyan
RGB(8,*)=[255,128,128]          ;Pink
RGB(9,*)=[128,128,255]          ;Violet
RGB(10,*)=[128,255,128]         ;Lime
RGB(11,*)=[255,128,0]           ;Brown
RGB(12,*)=[0,128,255]           ;Sea blue
RGB(13,*)=[255,0,128]           ;Red pink
RGB(14,*)=[128,0,255]           ;Purple
RGB(15,*)=[64,64,64]            ;Heave grey
RGB(16,*)=[128,128,128]         ;Mid grey
RGB(17,*)=[192,192,192]         ;Light grey

TVLCT, RGB

END
