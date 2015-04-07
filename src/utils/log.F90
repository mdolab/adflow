! File for open/closing an optional log file for SUmb
subroutine openLog(fileName)

  character*(*) :: fileName

  open(unit=6,file=filename,status="replace", access="sequential")

end subroutine openLog


subroutine closeLog()

  close(unit=6)

end subroutine closeLog
