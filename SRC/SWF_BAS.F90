!==================================================================================================
! 子程序包: SWF_BAS (Basic Package for Soil Water Flow Process)
! 主要功能: 计算土壤水分运动参数, 设置边界条件, 输出模拟结果等
! 创建日期: 2020年7月3日
! 修改日志:
!          2020年8月31日, 重构了 SWF_BAS 子程序包
!==================================================================================================


SUBROUTINE SWF_BAS_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: BAS, ONLY : U_LOG, U_IN_BAS, S_IN_BAS
  USE :: BAS, ONLY : APPNAM, PRJNAM, ISCENE, LCROP
  USE :: BAS, ONLY : NCEL, NPER
  USE :: SWF, ONLY : L_SWF_STO, L_SWF_CCF, L_SWF_GBB, L_SWF_GTB, L_SWF_INF, L_SWF_EVP, L_SWF_RWU
  USE :: SWF, ONLY : L_SWF_TMA
  USE :: SWF, ONLY : U_SWF_PFL, S_SWF_PFL, U_SWF_PER, S_SWF_PER, U_SWF_OBS, S_SWF_OBS
  USE :: SWF, ONLY : ISWFINI
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_BAS 文件读取数据
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('ISWFINI')) THEN
    CALL RDSINT('ISWFINI', ISWFINI)
  ELSE
    ISWFINI = 0
  END IF
  IF (RDINQR('S_SWF_PFL')) THEN
    CALL RDSCHA('S_SWF_PFL', S_SWF_PFL)
  ELSE
    S_SWF_PFL = 'SWF_PFL.OUT'
  END IF
  IF (RDINQR('S_SWF_PER')) THEN
    CALL RDSCHA('S_SWF_PER', S_SWF_PER)
  ELSE
    S_SWF_PER = 'SWF_PER.OUT'
  END IF
  IF (RDINQR('S_SWF_OBS')) THEN
    CALL RDSCHA('S_SWF_OBS', S_SWF_OBS)
  ELSE
    S_SWF_OBS = 'SWF_OBS.OUT'
  END IF

  ! 打开 S_SWF_PFL 文件, 输出标识信息
  OPEN(UNIT = U_SWF_PFL, FILE = S_SWF_PFL)
  WRITE(U_SWF_PFL, "(A)") APPNAM
  WRITE(U_SWF_PFL, "()")
  WRITE(U_SWF_PFL, "('土壤水分运动剖面数据输出文件')")
  WRITE(U_SWF_PFL, "()")
  WRITE(U_SWF_PFL, "('项目名称: ', A)") TRIM(PRJNAM)
  WRITE(U_SWF_PFL, "('应力期数量: ', I12)") NPER
  WRITE(U_SWF_PFL, "('单元体数量: ', I12)") NCEL
  WRITE(U_SWF_PFL, "()")

  ! 打开 S_SWF_PER 文件, 输出标识信息
  OPEN(UNIT = U_SWF_PER, FILE = S_SWF_PER)
  WRITE(U_SWF_PER, "(A)") APPNAM
  WRITE(U_SWF_PER, "()")
  WRITE(U_SWF_PER, "('土壤水分运动应力期数据输出文件')")
  WRITE(U_SWF_PER, "()")
  WRITE(U_SWF_PER, "('项目名称: ', A)") TRIM(PRJNAM)
  WRITE(U_SWF_PER, "('应力期数量: ', I12)") NPER
  WRITE(U_SWF_PER, "()")

  ! 打开 S_SWF_OBS 文件, 输出标识信息
  OPEN(UNIT = U_SWF_OBS, FILE = S_SWF_OBS)
  WRITE(U_SWF_OBS, '(A)') APPNAM
  WRITE(U_SWF_OBS, "()")
  WRITE(U_SWF_OBS, "('土壤水分运动观测点数据输出文件')")
  WRITE(U_SWF_OBS, "()")
  WRITE(U_SWF_OBS, "('项目名称: ', A)") TRIM(PRJNAM)
  WRITE(U_SWF_OBS, "('应力期数量: ', I12)") NPER
  WRITE(U_SWF_OBS, "()")

  L_SWF_STO = .TRUE.
  L_SWF_CCF = .TRUE.
  L_SWF_GBB = .TRUE.
  L_SWF_GTB = (ISCENE == 2)
  L_SWF_INF = (ISCENE == 1)
  L_SWF_EVP = (ISCENE == 1)
  L_SWF_RWU = LCROP
  L_SWF_TMA = .TRUE.

  RETURN
END SUBROUTINE SWF_BAS_DF


SUBROUTINE SWF_BAS_AL
!**************************************************************************************************
! 为数组分配内存并初始化
!**************************************************************************************************
  USE :: BAS, ONLY : NMAT, NCEL
  USE :: SWF, ONLY : ORES, OSAT, HYCONSAT, MVGA, MVGN, MVGL
  USE :: SWF, ONLY : HYCAP, HYCON, HCOF, RHS
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(ORES))     ALLOCATE(ORES(NMAT),     SOURCE = 0.0)
  IF (.NOT. ALLOCATED(OSAT))     ALLOCATE(OSAT(NMAT),     SOURCE = 0.0)
  IF (.NOT. ALLOCATED(HYCONSAT)) ALLOCATE(HYCONSAT(NMAT), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(MVGA))     ALLOCATE(MVGA(NMAT),     SOURCE = 0.0)
  IF (.NOT. ALLOCATED(MVGN))     ALLOCATE(MVGN(NMAT),     SOURCE = 0.0)
  IF (.NOT. ALLOCATED(MVGL))     ALLOCATE(MVGL(NMAT),     SOURCE = 0.0)

  IF (.NOT. ALLOCATED(HYCAP))    ALLOCATE(HYCAP(NCEL),    SOURCE = 0.0)
  IF (.NOT. ALLOCATED(HYCON))    ALLOCATE(HYCON(NCEL),    SOURCE = 0.0)
  IF (.NOT. ALLOCATED(HCOF))     ALLOCATE(HCOF(NCEL),     SOURCE = 0.0)
  IF (.NOT. ALLOCATED(RHS))      ALLOCATE(RHS(NCEL),      SOURCE = 0.0)

  RETURN
END SUBROUTINE SWF_BAS_AL


SUBROUTINE SWF_BAS_RP
!**************************************************************************************************
! 为数组读取数据并执行准备工作
!**************************************************************************************************
  USE :: BAS, ONLY : U_LOG, U_IN_SOL, S_IN_SOL
  USE :: BAS, ONLY : NMAT
  USE :: SWF, ONLY : ORES, OSAT, HYCONSAT, MVGA, MVGN, MVGL
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_SOL 文件读取数据
  CALL RDINIT(U_IN_SOL, U_LOG, S_IN_SOL)
  IF (RDINQR('ORES')) THEN
    CALL RDFREA('ORES', ORES, NMAT, NMAT)
  ELSE
    ORES(:) = 0.078
  END IF
  IF (RDINQR('OSAT')) THEN
    CALL RDFREA('OSAT', OSAT, NMAT, NMAT)
  ELSE
    OSAT(:) = 0.43
  END IF
  IF (RDINQR('HYCONSAT')) THEN
    CALL RDFREA('HYCONSAT', HYCONSAT, NMAT, NMAT)
  ELSE
    HYCONSAT(:) = 24.96
  END IF
  IF (RDINQR('MVGA')) THEN
    CALL RDFREA('MVGA', MVGA, NMAT, NMAT)
  ELSE
    MVGA(:) = 0.036
  END IF
  IF (RDINQR('MVGN')) THEN
    CALL RDFREA('MVGN', MVGN, NMAT, NMAT)
  ELSE
    MVGN(:) = 1.56
  END IF
  IF (RDINQR('MVGL')) THEN
    CALL RDFREA('MVGL', MVGL, NMAT, NMAT)
  ELSE
    MVGL(:) = 0.5
  END IF

  RETURN
END SUBROUTINE SWF_BAS_RP


SUBROUTINE SWF_BAS_FM
!**************************************************************************************************
! 计算土壤水分运动参数并重置 HCOF 和 RHS
!**************************************************************************************************
  USE :: BAS, ONLY : NCEL, IMAT
  USE :: SWF, ONLY : ORES, OSAT, HYCONSAT, MVGA, MVGN, MVGL
  USE :: SWF, ONLY : HYCAP, HYCON, HNEW, ONEW, HCOF, RHS
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I, J

  DO I = 1, NCEL
    J = IMAT(I)
    CALL U_O_MVG(ONEW(I), HNEW(I), ORES(J), OSAT(J), MVGA(J), MVGN(J))
    CALL U_HYCAP_MVG(HYCAP(I), HNEW(I), ORES(J), OSAT(J), MVGA(J), MVGN(J))
    CALL U_HYCON_MVG(HYCON(I), HNEW(I), HYCONSAT(J), MVGA(J), MVGN(J), MVGL(J))
    HCOF(I) = 0.0
    RHS(I) = 0.0
  END DO

  RETURN
END SUBROUTINE SWF_BAS_FM


SUBROUTINE SWF_BAS_OT
!**************************************************************************************************
! 输出土壤水分运动模拟结果
!**************************************************************************************************
  USE :: BAS, ONLY : NCEL, DELV, LOCV, IOBS
  USE :: BAS, ONLY : PERLEN, IPRNT, TOTIM, JPER
  USE :: SWF, ONLY : U_SWF_PFL, U_SWF_PER, U_SWF_OBS
  USE :: SWF, ONLY : HINI, HNEW, OINI, ONEW
  USE :: SWF, ONLY : WBOT, WTOP, WSURF, WRNOF, WINF, WEPOT, WEACT, WTPOT, WTACT
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I

  IF (JPER == 1) THEN

    ! 将标题输出到 U_SWF_PFL 文件
    WRITE(U_SWF_PFL, "('      PERIOD')", ADVANCE = 'NO')
    WRITE(U_SWF_PFL, "('        TIME')", ADVANCE = 'NO')
    WRITE(U_SWF_PFL, "('        CELL')", ADVANCE = 'NO')
    WRITE(U_SWF_PFL, "('       DEPTH')", ADVANCE = 'NO')
    IF (ALLOCATED(HNEW)) WRITE(U_SWF_PFL, "('        HEAD')", ADVANCE = 'NO')
    IF (ALLOCATED(ONEW)) WRITE(U_SWF_PFL, "('       MOIST')", ADVANCE = 'NO')
    WRITE(U_SWF_PFL, "()", ADVANCE = 'YES')

    ! 将标题输出到 U_SWF_PER 文件
    WRITE(U_SWF_PER, "('      PERIOD')", ADVANCE = 'NO')
    WRITE(U_SWF_PER, "('      PERLEN')", ADVANCE = 'NO')
    WRITE(U_SWF_PER, "('        TIME')", ADVANCE = 'NO')
    IF (ALLOCATED(ONEW))  WRITE(U_SWF_PER, "('        STOR')", ADVANCE = 'NO')
    IF (ALLOCATED(WBOT))  WRITE(U_SWF_PER, "('        WBOT')", ADVANCE = 'NO')
    IF (ALLOCATED(WTOP))  WRITE(U_SWF_PER, "('        WTOP')", ADVANCE = 'NO')
    IF (ALLOCATED(WSURF)) WRITE(U_SWF_PER, "('       WSURF')", ADVANCE = 'NO')
    IF (ALLOCATED(WRNOF)) WRITE(U_SWF_PER, "('       WRNOF')", ADVANCE = 'NO')
    IF (ALLOCATED(WINF))  WRITE(U_SWF_PER, "('        WINF')", ADVANCE = 'NO')
    IF (ALLOCATED(WEPOT)) WRITE(U_SWF_PER, "('       WEPOT')", ADVANCE = 'NO')
    IF (ALLOCATED(WEACT)) WRITE(U_SWF_PER, "('       WEACT')", ADVANCE = 'NO')
    IF (ALLOCATED(WTPOT)) WRITE(U_SWF_PER, "('       WTPOT')", ADVANCE = 'NO')
    IF (ALLOCATED(WTACT)) WRITE(U_SWF_PER, "('       WTACT')", ADVANCE = 'NO')

    IF (ALLOCATED(WBOT))  WRITE(U_SWF_PER, "('      C_WBOT')", ADVANCE = 'NO')
    IF (ALLOCATED(WTOP))  WRITE(U_SWF_PER, "('      C_WTOP')", ADVANCE = 'NO')
    IF (ALLOCATED(WSURF)) WRITE(U_SWF_PER, "('     C_WSURF')", ADVANCE = 'NO')
    IF (ALLOCATED(WRNOF)) WRITE(U_SWF_PER, "('     C_WRNOF')", ADVANCE = 'NO')
    IF (ALLOCATED(WINF))  WRITE(U_SWF_PER, "('      C_WINF')", ADVANCE = 'NO')
    IF (ALLOCATED(WEPOT)) WRITE(U_SWF_PER, "('     C_WEPOT')", ADVANCE = 'NO')
    IF (ALLOCATED(WEACT)) WRITE(U_SWF_PER, "('     C_WEACT')", ADVANCE = 'NO')
    IF (ALLOCATED(WTPOT)) WRITE(U_SWF_PER, "('     C_WTPOT')", ADVANCE = 'NO')
    IF (ALLOCATED(WTACT)) WRITE(U_SWF_PER, "('     C_WTACT')", ADVANCE = 'NO')
    WRITE(U_SWF_PER, "()", ADVANCE = 'YES')

    ! 将标题输出到 U_SWF_OBS 文件
    WRITE(U_SWF_OBS, "('        IOBS')", ADVANCE = 'NO')
    WRITE(U_SWF_OBS, "('        CELL')", ADVANCE = 'NO')
    WRITE(U_SWF_OBS, "('       DEPTH')", ADVANCE = 'NO')
    WRITE(U_SWF_OBS, "('      PERIOD')", ADVANCE = 'NO')
    WRITE(U_SWF_OBS, "('        TIME')", ADVANCE = 'NO')
    IF (ALLOCATED(HNEW)) WRITE(U_SWF_OBS, "('        HEAD')", ADVANCE = 'NO')
    IF (ALLOCATED(ONEW)) WRITE(U_SWF_OBS, "('       MOIST')", ADVANCE = 'NO')
    WRITE(U_SWF_OBS, "()", ADVANCE = 'YES')

    ! 将初始条件输出到 U_SWF_PFL 文件
    DO I = 1, NCEL
      WRITE(U_SWF_PFL, "(I12)",   ADVANCE = 'NO') 0
      WRITE(U_SWF_PFL, "(F12.4)", ADVANCE = 'NO') 0.0
      WRITE(U_SWF_PFL, "(I12)",   ADVANCE = 'NO') I
      IF (ALLOCATED(LOCV)) WRITE(U_SWF_PFL, "(F12.2)",  ADVANCE = 'NO') LOCV(I, 2)
      IF (ALLOCATED(HINI)) WRITE(U_SWF_PFL, "(ES12.3)", ADVANCE = 'NO') HINI(I)
      IF (ALLOCATED(OINI)) WRITE(U_SWF_PFL, "(F12.4)",  ADVANCE = 'NO') OINI(I)
      WRITE(U_SWF_PFL, "()", ADVANCE = 'YES')
    END DO

  END IF

  ! 将模拟结果输出到 U_SWF_PFL 文件
  IF ((IPRNT(JPER) == 1) .OR. (IPRNT(JPER) == 2) .OR. (IPRNT(JPER) == 3)) THEN
    DO I = 1, NCEL
      WRITE(U_SWF_PFL, "(I12)",   ADVANCE = 'NO') JPER
      WRITE(U_SWF_PFL, "(F12.4)", ADVANCE = 'NO') TOTIM
      WRITE(U_SWF_PFL, "(I12)",   ADVANCE = 'NO') I
      IF (ALLOCATED(LOCV)) WRITE(U_SWF_PFL, "(F12.2)",  ADVANCE = 'NO') LOCV(I, 2)
      IF (ALLOCATED(HNEW)) WRITE(U_SWF_PFL, "(ES12.3)", ADVANCE = 'NO') HNEW(I)
      IF (ALLOCATED(ONEW)) WRITE(U_SWF_PFL, "(F12.4)",  ADVANCE = 'NO') ONEW(I)
      WRITE(U_SWF_PFL, "()", ADVANCE = 'YES')
    END DO
  END IF

  ! 将模拟结果输出到 U_SWF_PER 文件
  IF ((IPRNT(JPER) == 1) .OR. (IPRNT(JPER) == 2) .OR. (IPRNT(JPER) == 4)) THEN
    WRITE(U_SWF_PER, "(I12)",   ADVANCE = 'NO') JPER
    WRITE(U_SWF_PER, "(F12.4)", ADVANCE = 'NO') PERLEN(JPER)
    WRITE(U_SWF_PER, "(F12.4)", ADVANCE = 'NO') TOTIM
    IF (ALLOCATED(ONEW))  WRITE(U_SWF_PER, "(F12.4)",  ADVANCE = 'NO') SUM(ONEW(:) * DELV(:))
    IF (ALLOCATED(WBOT))  WRITE(U_SWF_PER, "(ES12.3)", ADVANCE = 'NO') WBOT(JPER)
    IF (ALLOCATED(WTOP))  WRITE(U_SWF_PER, "(ES12.3)", ADVANCE = 'NO') WTOP(JPER)
    IF (ALLOCATED(WSURF)) WRITE(U_SWF_PER, "(ES12.3)", ADVANCE = 'NO') WSURF(JPER)
    IF (ALLOCATED(WRNOF)) WRITE(U_SWF_PER, "(ES12.3)", ADVANCE = 'NO') WRNOF(JPER)
    IF (ALLOCATED(WINF))  WRITE(U_SWF_PER, "(ES12.3)", ADVANCE = 'NO') WINF(JPER)
    IF (ALLOCATED(WEPOT)) WRITE(U_SWF_PER, "(ES12.3)", ADVANCE = 'NO') WEPOT(JPER)
    IF (ALLOCATED(WEACT)) WRITE(U_SWF_PER, "(ES12.3)", ADVANCE = 'NO') WEACT(JPER)
    IF (ALLOCATED(WTPOT)) WRITE(U_SWF_PER, "(ES12.3)", ADVANCE = 'NO') WTPOT(JPER)
    IF (ALLOCATED(WTACT)) WRITE(U_SWF_PER, "(ES12.3)", ADVANCE = 'NO') WTACT(JPER)

    IF (ALLOCATED(WBOT))  WRITE(U_SWF_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WBOT(1:JPER))
    IF (ALLOCATED(WTOP))  WRITE(U_SWF_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WTOP(1:JPER))
    IF (ALLOCATED(WSURF)) WRITE(U_SWF_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WSURF(1:JPER))
    IF (ALLOCATED(WRNOF)) WRITE(U_SWF_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WRNOF(1:JPER))
    IF (ALLOCATED(WINF))  WRITE(U_SWF_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WINF(1:JPER))
    IF (ALLOCATED(WEPOT)) WRITE(U_SWF_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WEPOT(1:JPER))
    IF (ALLOCATED(WEACT)) WRITE(U_SWF_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WEACT(1:JPER))
    IF (ALLOCATED(WTPOT)) WRITE(U_SWF_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WTPOT(1:JPER))
    IF (ALLOCATED(WTACT)) WRITE(U_SWF_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WTACT(1:JPER))
    WRITE(U_SWF_PER, "()", ADVANCE = 'YES')
  END IF

  ! 将模拟结果输出到 U_SWF_OBS 文件
  IF ((IPRNT(JPER) == 1) .OR. (IPRNT(JPER) == 3) .OR. (IPRNT(JPER) == 4)) THEN
    DO I = 1, NCEL
      IF (IOBS(I) > 0) THEN
        WRITE(U_SWF_OBS, "(I12)",   ADVANCE = 'NO') IOBS(I)
        WRITE(U_SWF_OBS, "(I12)",   ADVANCE = 'NO') I
        WRITE(U_SWF_OBS, "(F12.2)", ADVANCE = 'NO') LOCV(I, 2)
        WRITE(U_SWF_OBS, "(I12)",   ADVANCE = 'NO') JPER
        WRITE(U_SWF_OBS, "(F12.4)", ADVANCE = 'NO') TOTIM
        IF (ALLOCATED(HNEW)) WRITE(U_SWF_OBS, "(ES12.3)", ADVANCE = 'NO') HNEW(I)
        IF (ALLOCATED(ONEW)) WRITE(U_SWF_OBS, "(F12.4)",  ADVANCE = 'NO') ONEW(I)
        WRITE(U_SWF_OBS, "()", ADVANCE = 'YES')
      END IF
    END DO
  END IF

  RETURN
END SUBROUTINE SWF_BAS_OT


SUBROUTINE SWF_BAS_DA
!**************************************************************************************************
! 释放数组内存并关闭文件
!**************************************************************************************************
  USE :: SWF, ONLY : U_SWF_PFL, U_SWF_PER, U_SWF_OBS
  USE :: SWF, ONLY : ORES, OSAT, HYCONSAT, MVGA, MVGN, MVGL
  USE :: SWF, ONLY : HYCAP, HYCON, HCOF, RHS
  IMPLICIT NONE

  ! 释放数组变量内存
  IF (ALLOCATED(ORES))     DEALLOCATE(ORES)
  IF (ALLOCATED(OSAT))     DEALLOCATE(OSAT)
  IF (ALLOCATED(HYCONSAT)) DEALLOCATE(HYCONSAT)
  IF (ALLOCATED(MVGA))     DEALLOCATE(MVGA)
  IF (ALLOCATED(MVGN))     DEALLOCATE(MVGN)
  IF (ALLOCATED(MVGL))     DEALLOCATE(MVGL)

  IF (ALLOCATED(HYCAP))    DEALLOCATE(HYCAP)
  IF (ALLOCATED(HYCON))    DEALLOCATE(HYCON)
  IF (ALLOCATED(HCOF))     DEALLOCATE(HCOF)
  IF (ALLOCATED(RHS))      DEALLOCATE(RHS)

  ! 关闭输出文件
  CLOSE(U_SWF_PFL)
  CLOSE(U_SWF_PER)
  CLOSE(U_SWF_OBS)

  RETURN
END SUBROUTINE SWF_BAS_DA
