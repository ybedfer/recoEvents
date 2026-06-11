
ROOT_GENERATE_DICTIONARY(recoEventsDict recoEvents.h LINKDEF recoEventsLinkDef.h)

add_library(recoEvents SHARED recoEvents.C recoEventsDict.cxx geometry.C Geometry.h Geometry.cc )
target_link_libraries(recoEvents ${ROOT_LIBRARIES})
