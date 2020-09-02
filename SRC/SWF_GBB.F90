!==================================================================================================
! 子程序包: SWF_GBB (General Bottom Boundary Condition Package for Soil Water Flow Process)
! 主要功能: 模拟底部边界的水分运动
! 创建日期: 2020年7月3日
! 修改日志:
!          2020年8月31日, 重构了 SWF_GBB 子程序包
!==================================================================================================


SUBROUTINE SWF_GBB_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: BAS, ONLY : U_LOG, U_IN_BAS, S_IN_BAS
  USE :: SWF, ONLY : ISWFBOT, TLAGBOT
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_BAS 文件读取数据
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('ISWFBOT')) THEN
    CALL RDSINT('ISWFBOT', ISWFBOT)
  ELSE
    ISWFBOT = 0
  END IF
  IF (RDINQR('TLAGBOT')) THEN
    CALL RDSREA('TLAGBOT', TLAGBOT)
  ELSE
    TLAGBOT = 1.0
  END IF

  RETURN
END SUBROUTINE SWF_GBB_DF


SUBROUTINE SWF_GBB_AL
!**************************************************************************************************
! 为数组分配内存并初始化
!**************************************************************************************************
  USE :: BAS, ONLY : NPER
  USE :: SWF, ONLY : HBOT, OBOT, WBOT
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(HBOT)) ALLOCATE(HBOT(NPER), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(OBOT)) ALLOCATE(OBOT(NPER), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WBOT)) ALLOCATE(WBOT(NPER), SOURCE = 0.0)

  RETURN
END SUBROUTINE SWF_GBB_AL


SUBROUTINE SWF_GBB_RP
!**************************************************************************************************
! 为数组读取数据并执行准备工作
!**************************************************************************************************
  USE :: BAS, ONLY : U_LOG, U_IN_PER, S_IN_PER
  USE :: BAS, ONLY : NCEL, NPER, IMAT
  USE :: SWF, ONLY : ORES, OSAT, MVGA, MVGN
  USE :: SWF, ONLY : ISWFINI, HINI, OINI
  USE :: SWF, ONLY : ISWFBOT, HBOT, OBOT, WBOT
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 局部变量
  INTEGER :: I, J

  ! 从 S_IN_PER 文件读取数据
  CALL RDINIT(U_IN_PER, U_LOG, S_IN_PER)
  IF (RDINQR('HBOT')) THEN
    CALL RDFREA('HBOT', HBOT, NPER, NPER)
  ELSE
    HBOT(:) = HINI(NCEL)
  END IF
  IF (RDINQR('OBOT')) THEN
    CALL RDFREA('OBOT', OBOT, NPER, NPER)
  ELSE
    OBOT(:) = OINI(NCEL)
  END IF
  IF (RDINQR('WBOT')) THEN
    CALL RDFREA('WBOT', WBOT, NPER, NPER)
  ELSE
    WBOT(:) = 0.0
  END IF

  ! 根据初始条件类型换算水势和含水率
  J = IMAT(NCEL)
  DO I = 1, NPER
    IF (ISWFBOT == 1) THEN
      IF (ISWFINI == 0) THEN
        CALL U_O_MVG(OBOT(I), HBOT(I), ORES(J), OSAT(J), MVGA(J), MVGN(J))
      ELSE
        CALL U_H_MVG(HBOT(I), OBOT(I), ORES(J), OSAT(J), MVGA(J), MVGN(J))
      END IF
    END IF
  END DO

  RETURN
END SUBROUTINE SWF_GBB_RP


SUBROUTINE SWF_GBB_FM
!**************************************************************************************************
! 计算表示底部边界水分运动的方程组系数
!**************************************************************************************************
  USE :: BAS, ONLY : NCEL, IMAT, DELV, LOCV, DELT, JPER
  USE :: SWF, ONLY : ORES, OSAT, HYCONSAT, MVGA, MVGN, MVGL
  USE :: SWF, ONLY : ISWFINI, HYCON, HNEW, HCOF, RHS
  USE :: SWF, ONLY : ISWFBOT, TLAGBOT, HBOT, OBOT, WBOT
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I, J
  REAL :: DELV2
  REAL :: HYCONBOT, HYCON2

  I = NCEL
  J = IMAT(I)

  SELECT CASE (ISWFBOT)

  ! Dirichlet 条件
  CASE (1)
    IF (ISWFINI /= 0) THEN
      CALL U_H_MVG(HBOT(JPER), OBOT(JPER), ORES(J), OSAT(J), MVGA(J), MVGN(J))
    END IF
    CALL U_HYCON_MVG(HYCONBOT, HBOT(JPER), HYCONSAT(J), MVGA(J), MVGN(J), MVGL(J))
    HYCON2 = (HYCONBOT + HYCON(I)) * 0.5
    DELV2 = DELV(I) * 0.5
    HCOF(I) = HCOF(I) - DELT * HYCON2 / DELV2
    RHS(I) = RHS(I) - DELT * HYCON2 / DELV2 * HBOT(JPER) + DELT * HYCON2

  ! Neumann 条件
  CASE (2)
    RHS(I) = RHS(I) - DELT * WBOT(JPER)

  ! Cauchy 条件
  CASE (3)
    HCOF(I) = HCOF(I) - DELT / TLAGBOT
    RHS(I) = RHS(I) - DELT * (HBOT(JPER) - LOCV(I, 2)) / TLAGBOT

  ! 自由排水条件
  CASE (4)
    RHS(I) = RHS(I) + DELT * HYCON(I)

  ! 渗出面条件
  CASE (5)
    DELV2 = DELV(I) * 0.5
    HBOT(JPER) = HNEW(I) + DELV2
    IF (HBOT(JPER) >= 0.0) THEN
      HBOT(JPER) = 0.0
      CALL U_HYCON_MVG(HYCONBOT, HBOT(JPER), HYCONSAT(J), MVGA(J), MVGN(J), MVGL(J))
      HYCON2 = (HYCONBOT + HYCON(I)) * 0.5
      HCOF(I) = HCOF(I) - DELT * HYCON2 / DELV2
      RHS(I) = RHS(I) - DELT * (HYCON2 / DELV2 * HBOT(JPER) - HYCON2)
    END IF

  END SELECT

  RETURN
END SUBROUTINE SWF_GBB_FM


SUBROUTINE SWF_GBB_BD
!**************************************************************************************************
! 计算表示底部边界水分运动的均衡项
!**************************************************************************************************
  USE :: BAS, ONLY : NCEL, IMAT, DELV, LOCV, DELT, JPER
  USE :: SWF, ONLY : ORES, OSAT, HYCONSAT, MVGA, MVGN, MVGL
  USE :: SWF, ONLY : ISWFINI, HYCON, HNEW
  USE :: SWF, ONLY : ISWFBOT, TLAGBOT, HBOT, OBOT, WBOT, QBOT
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I, J
  REAL :: HYCONBOT, HYCON2
  REAL :: DELV2

  I = NCEL
  J = IMAT(I)


  SELECT CASE (ISWFBOT)

  ! Dirichlet 条件
  CASE (1)
    IF (ISWFINI /= 0) CALL U_H_MVG(HBOT(JPER), OBOT(JPER), ORES(J), OSAT(J), MVGA(J), MVGN(J))
    CALL U_HYCON_MVG(HYCONBOT, HBOT(JPER), HYCONSAT(J), MVGA(J), MVGN(J), MVGL(J))
    DELV2 = DELV(I) * 0.5
    HYCON2 = (HYCONBOT + HYCON(I)) * 0.5
    QBOT = DELT * HYCON2 * ((HBOT(JPER) - HNEW(I)) / DELV2 - 1.0)
    WBOT(JPER) = WBOT(JPER) + QBOT

  ! Neumann 条件
  CASE (2)
    QBOT = DELT * WBOT(JPER)

  ! Cauchy 条件
  CASE (3)
    QBOT = DELT * (HBOT(JPER) - LOCV(I, 2) - HNEW(I)) / TLAGBOT
    WBOT(JPER) = WBOT(JPER) + QBOT

  ! 自由排水条件
  CASE (4)
    QBOT = - DELT * HYCON(I)
    WBOT(JPER) = WBOT(JPER) + QBOT

  ! 渗出面条件
  CASE (5)
    IF (HBOT(JPER) >= 0.0) THEN
      CALL U_HYCON_MVG(HYCONBOT, HBOT(JPER), HYCONSAT(J), MVGA(J), MVGN(J), MVGL(J))
      DELV2 = DELV(I) * 0.5
      HYCON2 = (HYCONBOT + HYCON(I)) * 0.5
      QBOT = DELT * HYCON2 * ((HBOT(JPER) - HNEW(I)) / DELV2 - 1.0)
      WBOT(JPER) = WBOT(JPER) + QBOT
    END IF

  END SELECT

  RETURN
END SUBROUTINE SWF_GBB_BD


SUBROUTINE SWF_GBB_DA
!**************************************************************************************************
! 释放数组内存
!**************************************************************************************************
  USE :: SWF, ONLY : HBOT, OBOT, WBOT
  IMPLICIT NONE

  IF (ALLOCATED(HBOT)) DEALLOCATE(HBOT)
  IF (ALLOCATED(OBOT)) DEALLOCATE(OBOT)
  IF (ALLOCATED(WBOT)) DEALLOCATE(WBOT)

  RETURN
end subroutine SWF_GBB_DA
