#include <px4_platform_common/log.h>

__EXPORT int hello_world_main(int argc, char *argv[]);

int hello_world_main(int argc, char *argv[])
{
	PX4_INFO("Hello Word!");
	return OK;
}