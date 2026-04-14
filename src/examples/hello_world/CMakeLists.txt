px4_add_module(
	MODULE examples__hello_world
	MAIN hello_world
	STACK_MAIN 2000
	SRCS
		hello.c
	DEPENDS
)