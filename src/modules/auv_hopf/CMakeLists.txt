px4_add_module(
    MODULE modules__auv_hopf
    MAIN auv_hopf
    STACK_MAIN 6000
    SRCS
        auv_hopf_main.cpp
        auv_model.cpp
        hopf_observer.cpp
        wave_disturbance.cpp
        mpc_optim.cpp
        actuator_saturation.cpp
        telem_link.cpp
    INCLUDES
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PX4_SOURCE_DIR}/src/lib/grampc/include

    DEPENDS
        grampc
)
