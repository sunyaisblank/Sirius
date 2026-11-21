#include "packaged_software.h"

#include "sirius/backend/vulkan/vulkan_device.h"
#include "sirius/backend/vulkan/vulkan_portability.h"
#include "sirius/base/resource_locator.h"
#include "sirius/base/sha256.h"

#include "sealed_shaders.h"
#include <dlfcn.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <string_view>

using namespace sirius::backend;

namespace {

class VerifiedDriver {
  public:
    VerifiedDriver() = default;
    VerifiedDriver(const VerifiedDriver&) = delete;
    VerifiedDriver& operator=(const VerifiedDriver&) = delete;

    ~VerifiedDriver() {
        if (handle_ != nullptr) dlclose(handle_);
        if (file_ >= 0) close(file_);
    }

    std::expected<void, std::string> Open(const std::filesystem::path& path,
                                          std::string_view expected_sha256) {
        file_ = open(path.c_str(), O_RDONLY | O_CLOEXEC | O_NOFOLLOW);
        struct stat before {};
        if (file_ < 0 || fstat(file_, &before) != 0 || !S_ISREG(before.st_mode)) {
            return std::unexpected("cannot open driver as a regular file");
        }
        const std::string descriptor =
            "/proc/" + std::to_string(getpid()) + "/fd/" + std::to_string(file_);
        const auto digest = sirius::base::Sha256File(std::filesystem::path(descriptor));
        if (!digest || *digest != expected_sha256) {
            return std::unexpected("driver SHA-256 mismatch");
        }
        if (!Unchanged(before)) return std::unexpected("driver changed during verification");

        handle_ = dlopen(descriptor.c_str(), RTLD_NOW | RTLD_LOCAL);
        if (handle_ == nullptr) return std::unexpected(dlerror());
        void* symbol = dlsym(handle_, "vk_icdGetInstanceProcAddr");
        Dl_info owner{};
        struct stat mapped {};
        if (symbol == nullptr || dladdr(symbol, &owner) == 0 || owner.dli_fname == nullptr ||
            stat(owner.dli_fname, &mapped) != 0 || mapped.st_dev != before.st_dev ||
            mapped.st_ino != before.st_ino || !Unchanged(before)) {
            return std::unexpected("driver entry point is outside the verified file identity");
        }
        const void* llvm_symbol = dlsym(handle_, "LLVMGetVersion");
        Dl_info llvm_owner{};
        if (llvm_symbol == nullptr || dladdr(llvm_symbol, &llvm_owner) == 0 ||
            llvm_owner.dli_fname == nullptr) {
            return std::unexpected("cannot identify the loaded LLVM runtime");
        }
        std::error_code error;
        const auto llvm_path = std::filesystem::canonical(llvm_owner.dli_fname, error);
        if (error) return std::unexpected("cannot resolve the loaded LLVM runtime");
        const auto llvm_bytes = std::filesystem::file_size(llvm_path, error);
        if (error || llvm_bytes == 0)
            return std::unexpected("cannot determine the loaded LLVM runtime size");
        const auto llvm_digest = sirius::base::Sha256File(llvm_path);
        if (!llvm_digest || *llvm_digest != retained_program::kPackagedLlvmSha256) {
            return std::unexpected("LLVM runtime SHA-256 mismatch");
        }
        identity_ = PackagedVulkanProvider{
            {path.string(), static_cast<std::uint64_t>(before.st_size), *digest},
            {llvm_path.string(), llvm_bytes, *llvm_digest},
        };
        entry_ = reinterpret_cast<PFN_vkGetInstanceProcAddrLUNARG>(symbol);
        return {};
    }

    PFN_vkGetInstanceProcAddrLUNARG Entry() const { return entry_; }
    const std::optional<PackagedVulkanProvider>& Identity() const { return identity_; }

  private:
    bool Unchanged(const struct stat& before) const {
        struct stat now {};
        return fstat(file_, &now) == 0 && now.st_size == before.st_size &&
               now.st_mtim.tv_sec == before.st_mtim.tv_sec &&
               now.st_mtim.tv_nsec == before.st_mtim.tv_nsec &&
               now.st_ctim.tv_sec == before.st_ctim.tv_sec &&
               now.st_ctim.tv_nsec == before.st_ctim.tv_nsec;
    }

    int file_ = -1;
    void* handle_ = nullptr;
    PFN_vkGetInstanceProcAddrLUNARG entry_ = nullptr;
    std::optional<PackagedVulkanProvider> identity_;
};

struct Configuration {
    VerifiedDriver driver;
    std::optional<std::string> error;

    Configuration() {
        if (!detail::PackagedSoftwareRequested()) return;
        const auto path = sirius::base::ResolveResource("providers/linux-x86_64/libvulkan_lvp.so");
        if (!path) {
            error = "packaged software Vulkan driver is missing";
            return;
        }
        auto result = driver.Open(*path, retained_program::kPackagedDriverSha256);
        if (!result) error = std::move(result.error());
    }
};

const Configuration& BoundDriver() {
    static const Configuration configuration;
    return configuration;
}

}  // namespace

namespace sirius::backend {

std::optional<PackagedVulkanProvider> GetPackagedVulkanProvider() {
    return BoundDriver().driver.Identity();
}

}  // namespace sirius::backend

namespace sirius::backend::detail {

base::Expected<VkInstance> CreateProbedInstance(VkInstanceCreateInfo create_info) {
    const auto& configuration = BoundDriver();
    if (configuration.error)
        return base::Fail(base::ErrorDomain::kDevice, "open sealed Vulkan driver",
                          *configuration.error);
    const auto& driver = configuration.driver;
    if (driver.Entry() == nullptr) return CreateInstanceWithPortability(create_info);

    const VkDirectDriverLoadingInfoLUNARG verified_driver{
        .sType = VK_STRUCTURE_TYPE_DIRECT_DRIVER_LOADING_INFO_LUNARG,
        .pfnGetInstanceProcAddr = driver.Entry(),
    };
    const VkDirectDriverLoadingListLUNARG drivers{
        .sType = VK_STRUCTURE_TYPE_DIRECT_DRIVER_LOADING_LIST_LUNARG,
        .mode = VK_DIRECT_DRIVER_LOADING_MODE_EXCLUSIVE_LUNARG,
        .driverCount = 1,
        .pDrivers = &verified_driver,
    };
    const char* const enabled[] = {VK_LUNARG_DIRECT_DRIVER_LOADING_EXTENSION_NAME};
    create_info.pNext = &drivers;
    create_info.enabledExtensionCount = 1;
    create_info.ppEnabledExtensionNames = enabled;
    return CreateInstanceWithPortability(create_info);
}

bool UsesSealedRetainedModules(const DeviceInfo& info) {
    return BoundDriver().driver.Entry() != nullptr && info.kind == DeviceKind::kSoftware &&
           !RetainedUsesPortableArithmetic(info);
}

std::span<const std::uint32_t> SelectSealedRetainedModule(
    const DeviceInfo& info, std::span<const std::uint32_t> requested) {
    if (!UsesSealedRetainedModules(info)) return requested;
    using namespace retained_program;
    struct Module {
        std::string_view name, digest;
        std::span<const std::uint32_t> ordinary, sealed;
    };
    const Module known[] = {
        {"Transport", kTransportSealedSha256, kTransportShader, kTransportSealedShader},
        {"Endpoint", kEndpointSealedSha256, kEndpointShader, kEndpointSealedShader},
        {"DopriPhase", kDopriPhaseSealedSha256, kDopriPhaseShader, kDopriPhaseSealedShader},
    };
    for (const auto& module : known) {
        if (std::ranges::equal(requested, module.ordinary)) {
            std::cout << "module " << module.name << ' ' << module.sealed.size_bytes() << ' '
                      << module.digest << '\n';
            return module.sealed;
        }
    }
    return requested;
}

}  // namespace sirius::backend::detail
