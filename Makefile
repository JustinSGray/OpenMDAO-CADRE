#/bin/sh

Attitude = src/RungeKutta4.f90 \
	src/Attitude_Roll.f90 \
	src/Attitude_Attitude.f90 \
	src/Attitude_RotationMtx.f90 \
	src/Attitude_Angular.f90 \
	src/Attitude_Dynamics.f90 \
	src/Attitude_Torque.f90

Battery = src/RungeKutta4.f90 \
	src/Battery_SOC.f90

Comm = src/RungeKutta4.f90 \
	src/Comm_GSposEarth.f90 \
	src/Comm_EarthsSpin.f90 \
	src/Comm_AntRotation.f90 \
	src/Comm_Distance.f90 \
	src/Comm_LOS.f90 \
	src/Comm_BitRate.f90 \
	src/Comm_DataDownloaded.f90

Kinematics = src/Kinematics_RotationMtx.f90 \
	src/Kinematics_Rotate.f90 \
	src/Kinematics_Spherical.f90

KS = src/KSfunction.f90

Orbit = src/RungeKutta4.f90 \
	src/Orbit_Dynamics.f90

Power = src/Power_SolarPower.f90

RW = src/RungeKutta4.f90 \
	src/ReactionWheel_Dynamics.f90 \
	src/ReactionWheel_Motor.f90 \
	src/ReactionWheel_Power.f90

Sun = src/Sun_PositionECI.f90 \
	src/Sun_LOS.f90

Thermal = src/RungeKutta4.f90 \
	src/Thermal_Temperature.f90

default:	
	f2py --fcompiler=gnu95 -c ${Attitude} -m AttitudeLib
	f2py --fcompiler=gnu95 -c ${Battery} -m BatteryLib
	f2py --fcompiler=gnu95 -c ${Comm} -m CommLib
	f2py --fcompiler=gnu95 -c ${Kinematics} -m KinematicsLib
	f2py --fcompiler=gnu95 -c ${KS} -m KSLib
	f2py --fcompiler=gnu95 -c ${Orbit} -m OrbitLib
	f2py --fcompiler=gnu95 -c ${Power} -m PowerLib
	f2py --fcompiler=gnu95 -c ${RW} -m RWLib
	f2py --fcompiler=gnu95 -c ${Sun} -m SunLib
	f2py --fcompiler=gnu95 -c ${Thermal} -m ThermalLib
	-mv *.so ./CADRE/lib/
