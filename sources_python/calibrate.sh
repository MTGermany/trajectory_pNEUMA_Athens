#!/bin/bash

calibTraj 1 0 d8_0900_0930_road2_veh145 0
#exit 0
calibTraj 1 0 d8_0900_0930_road2_veh157 0
calibTraj 1 0 d8_0900_0930_road4_veh1085 0
calibTraj 1 0 d8_0900_0930_road4_veh1377 0

grep resulting d8*out
