      SUBROUTINE ADDREA (STRING,SIGLEN,R,FORM)
      IMPLICIT NONE

*     FORMAL_PARAMETERS:
      CHARACTER*(*) STRING,FORM
      INTEGER SIGLEN
      REAL R

**    local variables
      CHARACTER*30 TMPR
      CHARACTER*30 TMPFORM
      INTEGER IS,IL1,IL2,ISTART
      SAVE

      TMPFORM = '('//FORM//',A)'
      WRITE (TMPR,TMPFORM) R,'|'
      IL2 = INDEX (TMPR,'|')-1
      IS  = ISTART (TMPR)
      IL1 = LEN_TRIM (TMPR(1:IL2))

      IF (SIGLEN+(IL1-IS+1).GT.LEN (STRING)) CALL FATALERR
     &   ('ADDREA','internal error')

      STRING(SIGLEN+1:SIGLEN+IL1-IS+1) = TMPR(IS:IL1)
      SIGLEN = SIGLEN+IL1-IS+1

      RETURN
      END
