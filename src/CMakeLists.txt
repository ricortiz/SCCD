
set(SCCD_SOURCES
  DeformBVH.cpp
  orphan.cpp
  timing.cpp
  DeformModel.cpp
  ccd_vf.cpp
  intersect.cpp
  ccdAPI.cpp
  DeformBVH-ptr.cpp
)

add_library(sccd ${SCCD_SOURCES})
