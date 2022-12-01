#!/bin/bash

v0=13.785; T=1.522; s0=0.59; a=1.67; b=5.252
calibTraj 1 0 d8_0900_0930_road2_veh145 $v0 $T $s0 $a $b


v0=6.882;T=1.404;s0=0.421;a=0.969;b=4.701
calibTraj 1 0 d8_0900_0930_road2_veh157 $v0 $T $s0 $a $b

v0=9.047;T=2.160;s0=1.346;a=1.571;b=1.920
calibTraj 1 0 d8_0900_0930_road4_veh1377 $v0 $T $s0 $a $b


#calibTraj 1 0 d8_0900_0930_road4_veh1377
