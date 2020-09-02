!==================================================================================================
! 子程序包: BAS (Basic Package)
! 主要功能: 设置模拟场景, 运行模式, 模型规模, 时间离散, 空间离散等
! 创建日期: 2020年7月3日
! 修改日志:
!          2020年8月31日, 重构了 BAS 子程序包
!==================================================================================================


SUBROUTINE BAS_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: BAS, ONLY : U_LOG, S_LOG, U_IN_BAS, S_IN_BAS, S_IN_SOL, S_IN_PER
  USE :: BAS, ONLY : APPNAM, PRJNAM, ISCENE, LCROP, IRUNMOD
  USE :: BAS, ONLY : NMAT, NCEL, NPER
  USE :: BAS, ONLY : L_SWF_BAS
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 打开 S_LOG 文件, 输出标识信息
  OPEN(UNIT = U_LOG, FILE = S_LOG)
  WRITE(U_LOG, "(A)") APPNAM
  WRITE(U_LOG, "()")

  ! 从 S_IN_BAS 文件读取数据
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('PRJNAM')) THEN
    CALL RDSCHA('PRJNAM', PRJNAM)
  ELSE
    PRJNAM = 'This is a default project.'
  END IF
  IF (RDINQR('ISCENE')) THEN
    CALL RDSINT('ISCENE', ISCENE)
  ELSE
    ISCENE = 1
  END IF
  IF (RDINQR('LCROP')) THEN
    CALL RDSLOG('LCROP', LCROP)
  ELSE
    LCROP = .FALSE.
  END IF
  IF (RDINQR('IRUNMOD')) THEN
    CALL RDSINT('IRUNMOD', IRUNMOD)
  ELSE
    IRUNMOD = 1
  END IF
  IF (RDINQR('NMAT')) THEN
    CALL RDSINT('NMAT', NMAT)
  ELSE
    NMAT = 1
  END IF
  IF (RDINQR('NCEL')) THEN
    CALL RDSINT('NCEL', NCEL)
  ELSE
    NCEL = 100
  END IF
  IF (RDINQR('NPER')) THEN
    CALL RDSINT('NPER', NPER)
  ELSE
    NPER = 30
  END IF
  IF (RDINQR('S_IN_SOL')) THEN
    CALL RDSCHA('S_IN_SOL', S_IN_SOL)
  ELSE
    S_IN_SOL = 'SOIL.IN'
  END IF
  IF (RDINQR('S_IN_PER')) THEN
    CALL RDSCHA('S_IN_PER', S_IN_PER)
  ELSE
    S_IN_PER = 'PERIOD.IN'
  END IF

  L_SWF_BAS = .TRUE.

  RETURN
END SUBROUTINE BAS_DF


SUBROUTINE BAS_AL
!**************************************************************************************************
! 为数组分配内存并初始化
!**************************************************************************************************
  USE :: BAS, ONLY : NCEL, NPER
  USE :: BAS, ONLY : IMAT, IOBS, DELV, LOCV
  USE :: BAS, ONLY : PERLEN, NSTP, TSMLT, IPRNT
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(IMAT))   ALLOCATE(IMAT(NCEL),    SOURCE = 0)
  IF (.NOT. ALLOCATED(IOBS))   ALLOCATE(IOBS(NCEL),    SOURCE = 0)
  IF (.NOT. ALLOCATED(DELV))   ALLOCATE(DELV(NCEL),    SOURCE = 0.0)
  IF (.NOT. ALLOCATED(LOCV))   ALLOCATE(LOCV(NCEL, 3), SOURCE = 0.0)

  IF (.NOT. ALLOCATED(PERLEN)) ALLOCATE(PERLEN(NPER),  SOURCE = 0.0)
  IF (.NOT. ALLOCATED(NSTP))   ALLOCATE(NSTP(NPER),    SOURCE = 0)
  IF (.NOT. ALLOCATED(TSMLT))  ALLOCATE(TSMLT(NPER),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(IPRNT))  ALLOCATE(IPRNT(NPER),   SOURCE = 0)

  RETURN
END SUBROUTINE BAS_AL


SUBROUTINE BAS_RP
!**************************************************************************************************
! 为数组读取数据并执行准备工作
!**************************************************************************************************
  USE :: BAS, ONLY : U_LOG, U_IN_PER, S_IN_PER, U_IN_SOL, S_IN_SOL
  USE :: BAS, ONLY : IRUNMOD, IERRCOD, APPNAM, PRJNAM
  USE :: BAS, ONLY : NCEL, NPER
  USE :: BAS, ONLY : IMAT, DELV, IOBS, LOCV
  USE :: BAS, ONLY : PERLEN, NSTP, TSMLT, IPRNT, TOTIM
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 局部变量
  INTEGER :: I

  ! 从 S_IN_SOL 文件读取数据
  CALL RDINIT(U_IN_SOL, U_LOG, S_IN_SOL)
  IF (RDINQR('IMAT')) THEN
    CALL RDFINT('IMAT', IMAT, NCEL, NCEL)
  ELSE
    IMAT = 1
  END IF
  IF (RDINQR('DELV')) THEN
    CALL RDFREA('DELV', DELV, NCEL, NCEL)
  ELSE
    DELV(:) = 1.0
  END IF
  IF (RDINQR('IOBS')) THEN
    CALL RDFINT('IOBS', IOBS, NCEL, NCEL)
  ELSE
    IOBS(:) = 0
  END IF

  ! 从 S_IN_PER 文件读取数据
  CALL RDINIT(U_IN_PER, U_LOG, S_IN_PER)
  IF (RDINQR('PERLEN')) THEN
    CALL RDFREA('PERLEN', PERLEN, NPER, NPER)
  ELSE
    PERLEN(:) = 1.0
  END IF
  IF (RDINQR('NSTP')) THEN
    CALL RDFINT('NSTP', NSTP, NPER, NPER)
  ELSE
    NSTP(:) = 50
  END IF
  IF (RDINQR('TSMLT')) THEN
    CALL RDFREA('TSMLT', TSMLT, NPER, NPER)
  ELSE
    TSMLT(:) = 1.05
  END IF
  IF (RDINQR('IPRNT')) THEN
    CALL RDFINT('IPRNT', IPRNT, NPER, NPER)
  ELSE
    IPRNT(:) = 1
  END IF

  ! 在屏幕上显示标识信息
  IF (IRUNMOD == 1) THEN
    WRITE(*, "(A)") APPNAM
    WRITE(*, "()")
    WRITE(*, "('项目名称: ', A)") TRIM(PRJNAM)
    WRITE(*, "('应力期数量: ', I12)") NPER
    WRITE(*, "()")
  END IF

  ! 计算土壤单元体顶部, 中心和底部的位置
  LOCV(1, 1) = 0.0
  LOCV(1, 2) = DELV(1) * 0.5
  LOCV(1, 3) = DELV(1)
  DO I = 2, NCEL
    LOCV(I, 1) = LOCV(I - 1, 1) + DELV(I)
    LOCV(I, 2) = LOCV(I - 1, 2) + DELV(I)
    LOCV(I, 3) = LOCV(I - 1, 3) + DELV(I)
  END DO

  ! 初始化错误代码
  IERRCOD = 0

  ! 初始化模拟总时长
  TOTIM = 0.0

  RETURN
END SUBROUTINE BAS_RP


SUBROUTINE BAS_ST
!**************************************************************************************************
! 计算当前应力期的初始时段步长并输出当前应力期的信息
!**************************************************************************************************
  USE :: BAS, ONLY : U_LOG
  USE :: BAS, ONLY : IRUNMOD
  USE :: BAS, ONLY : PERLEN, NSTP, TSMLT, DELT, JPER
  IMPLICIT NONE

  ! 计算当前应力期的初始时段步长
  DELT = PERLEN(JPER) / FLOAT(NSTP(JPER))
  IF (TSMLT(JPER) /= 1.0) THEN
    DELT = PERLEN(JPER) * (1.0 - TSMLT(JPER)) / (1.0 - TSMLT(JPER) ** FLOAT(NSTP(JPER)))
  END IF

  ! 将当前应力期的序号输出到 S_LOG 文件, 并在屏幕上显示
  WRITE(U_LOG, "('  应力期: ', I12)") JPER
  IF (IRUNMOD == 1) WRITE(*, "('  应力期: ', I12)") JPER

  RETURN
END SUBROUTINE BAS_ST


SUBROUTINE BAS_AD
!**************************************************************************************************
! 计算当前时段步长和模拟总时长
!**************************************************************************************************
  USE :: BAS, ONLY : TSMLT, DELT, TOTIM, JPER, JSTP
  IMPLICIT NONE

  ! 计算当前时段步长
  IF (JSTP /= 1) DELT = DELT * TSMLT(JPER)

  ! 计算模拟总时长
  TOTIM = TOTIM + DELT

  RETURN
END SUBROUTINE BAS_AD


SUBROUTINE BAS_DA
!**************************************************************************************************
! 释放数组内存并关闭文件
!**************************************************************************************************
  USE :: BAS, ONLY : U_LOG, U_IN_BAS, U_IN_SOL, U_IN_PER
  USE :: BAS, ONLY : IMAT, DELV, IOBS, LOCV
  USE :: BAS, ONLY : PERLEN, NSTP, TSMLT, IPRNT
  IMPLICIT NONE

  ! 释放数组内存
  IF (ALLOCATED(IMAT))   DEALLOCATE(IMAT)
  IF (ALLOCATED(DELV))   DEALLOCATE(DELV)
  IF (ALLOCATED(IOBS))   DEALLOCATE(IOBS)
  IF (ALLOCATED(LOCV))   DEALLOCATE(LOCV)

  IF (ALLOCATED(PERLEN)) DEALLOCATE(PERLEN)
  IF (ALLOCATED(NSTP))   DEALLOCATE(NSTP)
  IF (ALLOCATED(TSMLT))  DEALLOCATE(TSMLT)
  IF (ALLOCATED(IPRNT))  DEALLOCATE(IPRNT)

 ! 删除临时文件
  CALL RDDTMP(U_IN_BAS)
  CALL RDDTMP(U_IN_SOL)
  CALL RDDTMP(U_IN_PER)

  ! 关闭输入文件
  CLOSE(U_IN_BAS)
  CLOSE(U_IN_SOL)
  CLOSE(U_IN_PER)

  ! 关闭日志文件
  CLOSE(U_LOG)

  RETURN
END SUBROUTINE BAS_DA
