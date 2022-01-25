set(SRCS
    rigidBody.cpp
    rigidBody_generalTensor.cpp
)
add_library(rigidBodyDynamics OBJECT ${SRCS})
target_link_libraries(rigidBodyDynamics PRIVATE quaternion)