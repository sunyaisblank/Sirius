#include "sirius/backend/vulkan/vulkan_device.h"

#include <array>
#include <fstream>
#include <iostream>
#include <span>
#include <stdexcept>
#include <vector>

using namespace sirius::backend;

template<class T>
T Check(sirius::base::Expected<T> value) {
    if (!value) throw std::runtime_error(value.error().Description());
    if constexpr (!std::is_void_v<T>) return std::move(*value);
}

template<class T>
std::vector<T> Read(const char* name) {
    std::ifstream file(name, std::ios::binary | std::ios::ate);
    const auto end = file.tellg();
    if (end <= 0 || end % sizeof(T) != 0)
        throw std::runtime_error("Wrong file size");
    std::vector<T> data(static_cast<std::size_t>(end)/sizeof(T));
    file.seekg(0);
    if (!file.read(reinterpret_cast<char*>(data.data()), end))
        throw std::runtime_error("Read input file");
    return data;
}

template<class T>
void Probe(const char* shader, const char* input, const char* output) {
    const auto inventory = Check(EnumerateVulkanDevices());
    auto device = Check(CreateVulkanDevice(Check(ResolveVulkanDeviceIndex(inventory))));
    std::cout << "device=" << device->Info().name << '\n';
    const auto kernel = Check(device->LoadKernel(Read<std::uint32_t>(shader)));
    const auto samples = Read<T>(input);
    if (samples.size() % 10 != 0 || samples.size() > 4096)
        throw std::runtime_error("Invalid sample roster");
    const auto count = static_cast<std::uint32_t>(samples.size()/10);
    std::vector<T> results(20*count);
    std::array<std::uint32_t,1> active{0};
    Check(device->SetBufferAllocationLimit(1024*1024));
    const std::array<BufferHandle,3> buffers{
        Check(device->CreateBuffer(sizeof(T)*samples.size(),BufferUsage::kStorage)),
        Check(device->CreateBuffer(sizeof(T)*results.size(),BufferUsage::kStorage)),
        Check(device->CreateBuffer(sizeof(active),BufferUsage::kStorage))
    };
    Check(device->WriteBuffer(buffers[0],std::as_bytes(std::span{samples})));
    Check(device->WriteBuffer(buffers[1],std::as_bytes(std::span{results})));
    Check(device->WriteBuffer(buffers[2],std::as_bytes(std::span{active})));
    DispatchTiming prepare;
    Check(device->Dispatch(kernel,buffers,1,1,1,&prepare));
    active[0]=count;
    Check(device->WriteBuffer(buffers[2],std::as_bytes(std::span{active})));
    DispatchTiming actual;
    Check(device->Dispatch(kernel,buffers,(count+7)/8,1,1,&actual));
    Check(device->ReadBuffer(buffers[1],std::as_writable_bytes(std::span{results})));
    std::cout << "actual_submit_ms=" << actual.submit_wait_ms << " samples=" << count << '\n';
    if (!(actual.submit_wait_ms > 0.0 && actual.submit_wait_ms <= 1000.0))
        throw std::runtime_error("Actual dispatch outside governor");
    if (!std::ofstream(output,std::ios::binary).write(
            reinterpret_cast<const char*>(results.data()),
            static_cast<std::streamsize>(sizeof(T)*results.size())))
        throw std::runtime_error("Write output");
}

int main(int argc, char** argv) try {
    if (argc != 5) throw std::runtime_error("Usage: chart-probe 32|64 shader input output");
    if(std::string_view(argv[1])=="32") Probe<float>(argv[2],argv[3],argv[4]);
    else if(std::string_view(argv[1])=="64") Probe<double>(argv[2],argv[3],argv[4]);
    else throw std::runtime_error("Expected 32 or 64");
    return 0;
} catch(const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
