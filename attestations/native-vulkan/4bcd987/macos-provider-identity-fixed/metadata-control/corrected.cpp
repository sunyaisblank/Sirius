#include <mvk_private_api.h>
#include <cstdio>
#include <cstring>
#define MVK_DIAGNOSTIC_CONTENT_ID "76cb8b212bc049261146b526d75ebddafe017b18fe1f4c2c7065c3024c4ffe9e"
int main() {
  VkPhysicalDeviceVulkan12Properties supportedProps12{};
  static_assert(sizeof(supportedProps12.driverInfo) == 256);
	std::snprintf(supportedProps12.driverInfo, sizeof(supportedProps12.driverInfo),
				  "%s+sirius-issue48-unforced.%s", MVK_VERSION_STRING, MVK_DIAGNOSTIC_CONTENT_ID);
  std::puts(supportedProps12.driverInfo);
  return std::strcmp(supportedProps12.driverInfo, "1.4.1+sirius-issue48-unforced.76cb8b212bc049261146b526d75ebddafe017b18fe1f4c2c7065c3024c4ffe9e") != 0;
}
