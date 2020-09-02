!==================================================================================================
! 子程序包: SWF_RWU (Root Water Uptake Pacakge for Soil Water Flow Process)
! 主要功能: 模拟作物根系吸水
! 创建日期: 2020年7月3日
! 修改日志:
!          2020年9月1日, 重构了 SWF_RWU 子程序包
!==================================================================================================


SUBROUTINE SWF_RWU_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: BAS, ONLY : U_LOG, U_IN_BAS, S_IN_BAS
  USE :: SWF, ONLY : IRWU, PRWUDIST, RWUCOMPCRT, RWUHCRT, RWUTCRT
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_BAS 文件读取数据
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('IRWU')) THEN
    CALL RDSINT('IRWU', IRWU)
  ELSE
    IRWU = 2
  END IF
  IF (RDINQR('PRWUDIST')) THEN
    CALL RDSREA('PRWUDIST', PRWUDIST)
  ELSE
    PRWUDIST = 5.0
  END IF
  IF (RDINQR('RWUCOMPCRT')) THEN
    CALL RDSREA('RWUCOMPCRT', RWUCOMPCRT)
  ELSE
    RWUCOMPCRT = 1.0
  END IF
  IF (RDINQR('RWUHCRT')) THEN
    CALL RDFREA('RWUHCRT', RWUHCRT, 5, 5)
  ELSE
    RWUHCRT(:) = (/0.0, -1.0, -300.0, -600.0, -15000.0/)
  END IF
  IF (RDINQR('RWUTCRT')) THEN
    CALL RDFREA('RWUTCRT', RWUTCRT, 2, 2)
  ELSE
    RWUTCRT(:) = (/0.5, 0.1/)
  END IF

  RETURN
END SUBROUTINE SWF_RWU_DF


SUBROUTINE SWF_RWU_AL
!**************************************************************************************************
! 为数组分配内存并初始化
!**************************************************************************************************
  USE :: BAS, ONLY : NPER, NCEL
  USE :: SWF, ONLY : RDEPTH, WTPOT, WTACT, RWUDIST, WRWUACT, QRWUACT
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(RDEPTH))  ALLOCATE(RDEPTH(NPER),  SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WTPOT))   ALLOCATE(WTPOT(NPER),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WTACT))   ALLOCATE(WTACT(NPER),   SOURCE = 0.0)

  IF (.NOT. ALLOCATED(RWUDIST)) ALLOCATE(RWUDIST(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WRWUACT)) ALLOCATE(WRWUACT(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(QRWUACT)) ALLOCATE(QRWUACT(NCEL), SOURCE = 0.0)

  RETURN
END SUBROUTINE SWF_RWU_AL


SUBROUTINE SWF_RWU_RP
!**************************************************************************************************
! 为数组读取数据并执行准备工作
!**************************************************************************************************
  USE :: BAS, ONLY : U_LOG, U_IN_PER, S_IN_PER, U_IN_SOL, S_IN_SOL
  USE :: BAS, ONLY : NPER, NCEL, LOCV
  USE :: SWF, ONLY : IRWU, WTPOT, RDEPTH, WRWUACT, RWUDIST, RWUCOMP
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 局部变量
  INTEGER :: I
  REAL :: S

  ! 从 S_IN_PER 文件读取数据
  CALL RDINIT(U_IN_PER, U_LOG, S_IN_PER)
  IF (RDINQR('WTPOT')) THEN
    CALL RDFREA('WTPOT', WTPOT, NPER, NPER)
  ELSE
    WTPOT(:) = 0.4
  END IF
  IF (RDINQR('RDEPTH')) THEN
    CALL RDFREA('RDEPTH', RDEPTH, NPER, NPER)
  ELSE
    RDEPTH(:) = LOCV(NCEL, 2) * 0.5
  END IF

  ! 从 S_IN_SOL 文件读取数据
  CALL RDINIT(U_IN_SOL, U_LOG, S_IN_SOL)
  IF (RDINQR('WRWUACT')) THEN
    CALL RDFREA('WRWUACT', WRWUACT, NCEL, NCEL)
  END IF
  IF (RDINQR('RWUDIST')) THEN
    CALL RDFREA('RWUDIST', RWUDIST, NCEL, NCEL)
  END IF

  ! 初始化根系吸水补偿系数
  RWUCOMP = 1.0

  ! 对输入的根系吸水分布进行归一化处理
  IF (IRWU == 1) THEN
    S = SUM(RWUDIST(:))
    DO I = 1, NCEL
      RWUDIST(I) = RWUDIST(I) / S
    END DO
  END IF

  RETURN
END SUBROUTINE SWF_RWU_RP


SUBROUTINE SWF_RWU_ST
!**************************************************************************************************
! 更新随应力期变化的数据
!**************************************************************************************************
  USE :: BAS, ONLY : NCEL, JPER, LOCV
  USE :: SWF, ONLY : IRWU, RDEPTH, RWUDIST, PRWUDIST
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I
  REAL :: R1, R2, SUM1

  ! 计算根系吸水分布
  IF (IRWU == 2) THEN
    DO I = 1, NCEL
      RWUDIST(I) = 0.0
      IF (RDEPTH(JPER) >= LOCV(I, 2)) THEN
        R1 = EXP(- PRWUDIST * LOCV(I, 1) / RDEPTH(JPER)) / (1.0 - EXP(- PRWUDIST))
        R2 = EXP(- PRWUDIST * LOCV(I, 3) / RDEPTH(JPER)) / (1.0 - EXP(- PRWUDIST))
        RWUDIST(I) = R1 - R2
      END IF
    END DO
  END IF

  RETURN
END SUBROUTINE SWF_RWU_ST


SUBROUTINE SWF_RWU_FM
!**************************************************************************************************
! 计算表示根系吸水的方程组系数
!**************************************************************************************************
  USE :: BAS, ONLY : NCEL, PERLEN, DELT, JPER
  USE :: SWF, ONLY : HNEW, RHS
  USE :: SWF, ONLY : IRWU, WTPOT, RWUDIST, RWUCOMP, RWUHCRT, RWUTCRT, WRWUACT, QRWUACT
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I
  REAL :: WSR      ! 根系吸水的水分胁迫系数

  IF (IRWU == 0) THEN
    RHS(I) = RHS(I) + DELT * WRWUACT(I)
  END IF

  IF ((IRWU == 1) .OR. (IRWU == 2)) THEN
    QRWUACT(:) = 0.0
    DO I = 1, NCEL
      CALL U_RWUWSR_FED(WSR, HNEW(I), WTPOT(JPER) / PERLEN(JPER), RWUHCRT, RWUTCRT)
      QRWUACT(I) = DELT * WTPOT(JPER) * RWUDIST(I) * WSR / RWUCOMP
      RHS(I) = RHS(I) + QRWUACT(I)
    END DO
  END IF

  RETURN
END SUBROUTINE SWF_RWU_FM


SUBROUTINE SWF_RWU_BD
!**************************************************************************************************
! 计算表示根系吸水的均衡项
!**************************************************************************************************
  USE :: BAS, ONLY : JPER, DELT
  USE :: SWF, ONLY : IRWU, WTPOT, WTACT, RWUDIST, RWUCOMPCRT, RWUCOMP, QRWUACT, WRWUACT
  IMPLICIT NONE

  IF ((IRWU == 1) .OR. (IRWU == 2)) THEN
    WRWUACT(:) = WRWUACT(:) + QRWUACT(:)
    WTACT(JPER) = WTACT(JPER) + SUM(QRWUACT(:))

    ! 更新根系吸水补偿系数
    RWUCOMP = 1.0
    IF (SUM(QRWUACT(:)) < DELT * WTPOT(JPER)) THEN
      IF (SUM(QRWUACT(:)) / (DELT * WTPOT(JPER)) > RWUCOMPCRT) THEN
        RWUCOMP = SUM(QRWUACT(:)) / SUM(DELT * WTPOT(JPER) * RWUDIST(:))
      END IF
    END IF
  END IF

  RETURN
END SUBROUTINE SWF_RWU_BD


SUBROUTINE SWF_RWU_DA
!**************************************************************************************************
! 释放数组内存
!**************************************************************************************************
  USE :: SWF, ONLY : RDEPTH, WTPOT, WTACT, RWUDIST, WRWUACT, QRWUACT
  IMPLICIT NONE

  IF (ALLOCATED(RDEPTH))  DEALLOCATE(RDEPTH)
  IF (ALLOCATED(WTPOT))   DEALLOCATE(WTPOT)
  IF (ALLOCATED(WTACT))   DEALLOCATE(WTACT)

  IF (ALLOCATED(RWUDIST)) DEALLOCATE(RWUDIST)
  IF (ALLOCATED(WRWUACT)) DEALLOCATE(WRWUACT)
  IF (ALLOCATED(QRWUACT)) DEALLOCATE(QRWUACT)

  RETURN
END SUBROUTINE SWF_RWU_DA

