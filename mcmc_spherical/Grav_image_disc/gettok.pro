function gettok,st,char
;+
; NAME:
;	GETTOK                                    
; PURPOSE:
;	Retrieve the first part of the string up to a specified character
; EXPLANATION:
;	GET TOKen - Retrieve first part of string until the character char 
;	is encountered.
;
; CALLING SEQUENCE:
;	token = gettok( st, char )
;
; INPUT:
;	char - character separating tokens, scalar string
;
; INPUT-OUTPUT:
;	st - (scalar) string to get token from (on output token is removed)
;
; OUTPUT:
;	token - scalar string value is returned 
;
; EXAMPLE:
;	If ST is 'abc=999' then gettok(ST,'=') would return
;	'abc' and ST would be left as '999' 
;
; NOTES:
;       A version of GETTOK that accepts vector strings is available for users 
;       of IDL V5.3 or later from  http://idlastro.gsfc.nasa.gov/ftp/v53/
; HISTORY
;	version 1  by D. Lindler APR,86
;	Remove leading blanks    W. Landsman (from JKF)    Aug. 1991
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;----------------------------------------------------------------------
  On_error,2                           ;Return to caller

; if char is a blank treat tabs as blanks

  tab = string(9b)
  while strpos(st,tab) GE 0 do begin    ;Search for tabs
	pos = strpos(st,tab)
	strput,st,' ',pos
  endwhile

  st = strtrim(st,1)              ;Remove leading blanks

; find character in string

  pos = strpos(st,char)
  if pos EQ -1 then begin         ;char not found?
	token = st
 	st = ''
 	return, token
  endif

; extract token

 token = strmid(st,0,pos)
 len = strlen(st)
 if pos EQ (len-1) then st = '' else st = strmid(st,pos+1,len-pos-1)

;  Return the result.

 return,token
 end

function numlines,file
;+
; NAME:
;     NUMLINES() 
; PURPOSE:
;     Return the number of lines in a file
;     This procedures became mostly obsolete in V5.6 with the introduction of
;     the FILE_LINES() procedure
; CALLING SEQUENCE:
;     nl = NUMLINES( filename )
; INPUT:
;     filename = name of file, scalar string
; OUTPUT:
;     nl = number of lines in the file, scalar longword
;          Set to -1 if the number of lines could not be determined
; METHOD:
;     If Unix then spawn to wc; otherwise read 1 line at a time and count
;
; PROCEDURE CALLS:
;     EXPAND_TILDE(), SPEC_DIR()
; MODIFICATION HISTORY:
;     W. Landsman                              February 1996
;     Use /bin/sh shell with wc under Unix     March 1997
;     Use EXPAND_TILDE() under Unix         September 1997
;     Converted to IDL V5.0   W. Landsman   September 1997
;-
 On_error,2

 if N_params() EQ 0 then begin
        print,'Syntax - nl = NUMLINES( file)'
        return,-1
 endif

 if !VERSION.RELEASE GE '5.6' then return,file_lines(file)
  nl = -1L
 openr,lun,file,/get_lun, ERROR = err
 if err NE 0 then begin
        if !VERSION.OS eq "vms" then file = spec_dir(file,'DAT') else $
        file = spec_dir(file)
        message,'ERROR - Unable to open file '+ file,/CON
        return,-1
 endif

 if !VERSION.OS_FAMILY EQ 'unix' then begin
         free_lun,lun
         if strpos(file,'~') GE 0 then file = expand_tilde(file)
         spawn,'wc -l < '+file, result, /sh    
         return,long(result[0])
 endif else begin                 ;=====>> Loop through file counting lines  
        On_ioerror,NOASCII
        nl = 0l
        tmp = ' '
         while not eof(lun) do begin
           readf,lun,tmp
           nl = nl + 1
         endwhile
         free_lun,lun
         return,nl
 endelse

NOASCII:
  message,'Error reading file ' + string(file),/CON
  return,-1
 end
