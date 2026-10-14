// Focused new-reconstruction comparator. No renderer or physics code.
#define STBI_ONLY_PNG
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>
#define TINYEXR_USE_MINIZ 1
#define TINYEXR_IMPLEMENTATION
#include <tinyexr.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using Bytes = std::vector<unsigned char>;
static_assert(sizeof(float) == 4 && std::numeric_limits<float>::is_iec559,
              "EXR domain requires IEEE binary32 float words");

void Require(bool ok, const std::string& why) {
    if (!ok) throw std::runtime_error(why);
}

std::string Quoted(const std::string& value) {
    std::ostringstream out;
    out << '"';
    for (unsigned char c : value) {
        if (c == '"' || c == '\\') out << '\\' << c;
        else if (c >= 32 && c < 127) out << c;
        else out << "?";
    }
    return out.str() + '"';
}

Bytes Read(const fs::path& path) {
    Require(fs::is_regular_file(path) && !fs::is_symlink(path), "not an ordinary input file");
    const auto size = fs::file_size(path);
    Require(size > 0 && size <= 256u * 1024u * 1024u, "input size outside bounded decoder domain");
    Bytes data(static_cast<size_t>(size));
    std::ifstream file(path, std::ios::binary);
    Require(bool(file.read(reinterpret_cast<char*>(data.data()), data.size())), "short input read");
    Require(file.peek() == std::char_traits<char>::eof(), "input grew during read");
    return data;
}

void Write(const fs::path& path, const Bytes& data) {
    Require(!fs::exists(path), "refuse to overwrite a fixture");
    std::ofstream file(path, std::ios::binary);
    Require(bool(file.write(reinterpret_cast<const char*>(data.data()), data.size())), "fixture write failed");
}

uint32_t Big32(const unsigned char* data) {
    return (uint32_t(data[0]) << 24) | (uint32_t(data[1]) << 16) |
           (uint32_t(data[2]) << 8) | uint32_t(data[3]);
}

void AppendBig32(Bytes& out, uint32_t value) {
    for (int shift : {24, 16, 8, 0}) out.push_back(static_cast<unsigned char>(value >> shift));
}

uint32_t Crc32(const unsigned char* data, size_t size) {
    uint32_t crc = ~uint32_t(0);
    for (size_t i = 0; i < size; ++i) {
        crc ^= data[i];
        for (int bit = 0; bit < 8; ++bit) crc = (crc >> 1) ^ (0xedb88320u & (0u - (crc & 1u)));
    }
    return ~crc;
}

struct Image {
    int width = 0, height = 0, channels = 0;
    std::string domain;
    std::vector<uint32_t> words;
};

size_t Samples(int width, int height, int channels) {
    Require(width > 0 && height > 0 && width <= 16384 && height <= 16384,
            "invalid expected dimensions");
    const uint64_t count = uint64_t(width) * uint64_t(height) * uint64_t(channels);
    Require(count <= 64u * 1024u * 1024u, "sample allocation exceeds bounded decoder domain");
    return static_cast<size_t>(count);
}

Image Png(const fs::path& path, int width, int height) {
    const auto file = Read(path);
    const unsigned char signature[8] = {137, 80, 78, 71, 13, 10, 26, 10};
    Require(file.size() >= 33 && std::equal(signature, signature + 8, file.begin()), "not PNG");
    Require(Big32(file.data() + 8) == 13 && std::memcmp(file.data() + 12, "IHDR", 4) == 0,
            "PNG must start with IHDR");
    Require(Big32(file.data() + 16) == uint32_t(width) && Big32(file.data() + 20) == uint32_t(height),
            "PNG dimensions differ from declared domain");
    Require(file[24] == 8 && file[25] == 6 && file[26] == 0 && file[27] == 0 && file[28] == 0,
            "PNG domain is noninterlaced RGBA8");
    bool ended = false, data_seen = false;
    size_t offset = 8;
    while (offset < file.size()) {
        Require(file.size() - offset >= 12, "truncated PNG chunk");
        const uint32_t length = Big32(file.data() + offset);
        Require(length <= file.size() - offset - 12, "PNG chunk exceeds file");
        const auto* type = file.data() + offset + 4;
        Require(Crc32(type, size_t(length) + 4) == Big32(file.data() + offset + length + 8), "PNG CRC differs");
        if (std::memcmp(type, "IDAT", 4) == 0) data_seen = true;
        if (std::memcmp(type, "IEND", 4) == 0) {
            Require(length == 0 && data_seen && offset + 12 == file.size(), "invalid PNG end");
            ended = true;
        }
        offset += size_t(length) + 12;
    }
    Require(ended, "missing PNG end");
    int actual_width = 0, actual_height = 0, channels = 0;
    std::unique_ptr<stbi_uc, decltype(&stbi_image_free)> pixels(
        stbi_load(path.string().c_str(), &actual_width, &actual_height, &channels, 4), stbi_image_free);
    Require(bool(pixels) && actual_width == width && actual_height == height && channels == 4,
            "stb PNG decode/dimensions/channels failed");
    Image image{width, height, 4, "RGBA8", {}};
    const size_t count = Samples(width, height, 4);
    image.words.reserve(count);
    for (size_t i = 0; i < count; ++i) image.words.push_back(pixels.get()[i]);
    return image;
}

struct Header {
    EXRHeader value;
    Header() { InitEXRHeader(&value); }
    ~Header() { FreeEXRHeader(&value); }
};

void ExrSuccess(int code, const char* error, const char* phase) {
    const std::string detail = error ? error : "";
    if (error) FreeEXRErrorMessage(error);
    Require(code == TINYEXR_SUCCESS, std::string(phase) + ": " + detail);
}

Image Exr(const fs::path& path, int width, int height) {
    (void)Read(path);
    EXRVersion version{};
    Require(ParseEXRVersionFromFile(&version, path.string().c_str()) == TINYEXR_SUCCESS,
            "EXR version failed");
    Require(!version.tiled && !version.multipart && !version.non_image, "EXR must be single scanline image");
    Header header;
    const char* error = nullptr;
    const int header_result = ParseEXRHeaderFromFile(&header.value, &version, path.string().c_str(), &error);
    ExrSuccess(header_result, error, "EXR header");
    const auto& h = header.value;
    Require(h.num_channels == 4 && h.compression_type == TINYEXR_COMPRESSIONTYPE_ZIP,
            "EXR domain requires four ZIP HALF channels");
    const std::array<std::string, 4> expected{"A", "B", "G", "R"};
    for (int i = 0; i < 4; ++i) {
        Require(h.channels[i].name == expected[i] && h.pixel_types[i] == TINYEXR_PIXELTYPE_HALF &&
                h.channels[i].x_sampling == 1 && h.channels[i].y_sampling == 1,
                "EXR channel name/type/sampling differs");
    }
    Require(h.data_window.min_x == 0 && h.data_window.min_y == 0 &&
            h.data_window.max_x == width - 1 && h.data_window.max_y == height - 1 &&
            h.display_window.min_x == 0 && h.display_window.min_y == 0 &&
            h.display_window.max_x == width - 1 && h.display_window.max_y == height - 1,
            "EXR windows differ from declared domain");
    float* raw = nullptr;
    int actual_width = 0, actual_height = 0;
    error = nullptr;
    const int code = LoadEXR(&raw, &actual_width, &actual_height, path.string().c_str(), &error);
    std::unique_ptr<float, decltype(&std::free)> pixels(raw, std::free);
    ExrSuccess(code, error, "EXR pixels");
    Require(bool(pixels) && actual_width == width && actual_height == height, "EXR dimensions differ");
    Image image{width, height, 4, "decoded-RGBA-IEEE754-binary32-words-from-HALF", {}};
    const size_t count = Samples(width, height, 4);
    image.words.resize(count);
    for (size_t i = 0; i < count; ++i) std::memcpy(&image.words[i], pixels.get() + i, 4);
    // No finite-value substitution, rounding, epsilon, or float arithmetic.
    return image;
}

Image Ppm(const fs::path& path, int width, int height) {
    const auto file = Read(path);
    size_t offset = 0;
    const auto token = [&]() {
        while (offset < file.size()) {
            if (std::isspace(file[offset])) { ++offset; continue; }
            if (file[offset] == '#') {
                while (offset < file.size() && file[offset] != '\n') ++offset;
                continue;
            }
            break;
        }
        const size_t first = offset;
        while (offset < file.size() && !std::isspace(file[offset]) && file[offset] != '#') ++offset;
        Require(offset > first, "missing PPM token");
        return std::string(file.begin() + first, file.begin() + offset);
    };
    Require(token() == "P6" && token() == std::to_string(width) && token() == std::to_string(height) &&
            token() == "255", "PPM header differs from declared P6 RGB8 domain");
    Require(offset < file.size() && std::isspace(file[offset]), "missing PPM payload delimiter");
    if (file[offset++] == '\r' && offset < file.size() && file[offset] == '\n') ++offset;
    const size_t count = Samples(width, height, 3);
    Require(file.size() - offset == count, "PPM payload is short or has trailing bytes");
    Image image{width, height, 3, "RGB8-P6", {}};
    image.words.reserve(count);
    for (; offset < file.size(); ++offset) image.words.push_back(file[offset]);
    return image;
}

Image Decode(const std::string& format, const fs::path& path, int width, int height) {
    (void)Samples(width, height, format == "ppm" ? 3 : 4);
    if (format == "png") return Png(path, width, height);
    if (format == "exr") return Exr(path, width, height);
    if (format == "ppm") return Ppm(path, width, height);
    throw std::runtime_error("unsupported format");
}

int Compare(const std::string& format, int width, int height, const fs::path& left, const fs::path& right) {
    const auto a = Decode(format, left, width, height);
    const auto b = Decode(format, right, width, height);
    Require(a.channels == b.channels && a.domain == b.domain && a.words.size() == b.words.size(),
            "decoded domains differ");
    size_t changed_pixels = 0, changed_samples = 0;
    for (size_t pixel = 0; pixel < size_t(width) * height; ++pixel) {
        bool changed = false;
        for (int channel = 0; channel < a.channels; ++channel) {
            const size_t index = pixel * a.channels + channel;
            if (a.words[index] != b.words[index]) { changed = true; ++changed_samples; }
        }
        if (changed) ++changed_pixels;
    }
    std::cout << "{\"equal\":" << (changed_pixels ? "false" : "true")
              << ",\"format\":" << Quoted(format) << ",\"sample_domain\":" << Quoted(a.domain)
              << ",\"width\":" << width << ",\"height\":" << height << ",\"channels\":" << a.channels
              << ",\"compared_pixels\":" << size_t(width) * height << ",\"compared_samples\":" << a.words.size()
              << ",\"mismatched_pixels\":" << changed_pixels << ",\"mismatched_samples\":" << changed_samples
              << ",\"encoded_files_equal\":" << (Read(left) == Read(right) ? "true" : "false") << "}\n";
    return changed_pixels ? 1 : 0;
}

void ExrFixture(const fs::path& path, bool metadata, bool perturb) {
    std::array<std::vector<float>, 4> planes;
    for (int channel = 0; channel < 4; ++channel) planes[channel].resize(6);
    for (size_t pixel = 0; pixel < 6; ++pixel) {
        planes[0][pixel] = 1.0f;                  // A
        planes[1][pixel] = 1.5f;                  // B
        planes[2][pixel] = 0.5f;                  // G
        planes[3][pixel] = 0.25f + pixel * 0.25f; // R
    }
    if (perturb) planes[2][0] = 0.50048828125f; // one represented HALF ULP
    EXRImage image;
    InitEXRImage(&image);
    image.width = 3; image.height = 2; image.num_channels = 4;
    unsigned char* channels[4];
    for (int i = 0; i < 4; ++i) channels[i] = reinterpret_cast<unsigned char*>(planes[i].data());
    image.images = channels;
    EXRHeader header;
    InitEXRHeader(&header);
    EXRChannelInfo info[4]{};
    const char* names[] = {"A", "B", "G", "R"};
    int types[4], requested[4];
    for (int i = 0; i < 4; ++i) {
        std::strcpy(info[i].name, names[i]); types[i] = TINYEXR_PIXELTYPE_FLOAT;
        info[i].x_sampling = info[i].y_sampling = 1;
        requested[i] = TINYEXR_PIXELTYPE_HALF;
    }
    header.channels = info; header.num_channels = 4;
    header.pixel_types = types; header.requested_pixel_types = requested;
    header.compression_type = TINYEXR_COMPRESSIONTYPE_ZIP;
    EXRAttribute attribute{};
    std::string note = "same pixels, extra metadata";
    if (metadata) {
        std::strcpy(attribute.name, "forensic_fixture_note"); std::strcpy(attribute.type, "string");
        attribute.value = reinterpret_cast<unsigned char*>(note.data()); attribute.size = int(note.size());
        header.custom_attributes = &attribute; header.num_custom_attributes = 1;
    }
    const char* error = nullptr;
    const int result = SaveEXRImageToFile(&image, &header, path.string().c_str(), &error);
    ExrSuccess(result, error, "EXR fixture write");
}

void Fixtures(const fs::path& folder) {
    Require(!fs::exists(folder) && fs::create_directory(folder), "fixture directory already exists or parent missing");
    Bytes rgba(24), rgb(18);
    for (size_t i = 0; i < rgba.size(); ++i) rgba[i] = i % 4 == 3 ? 255 : static_cast<unsigned char>((i * 37 + 11) % 250);
    for (size_t i = 0; i < rgb.size(); ++i) rgb[i] = static_cast<unsigned char>((i * 37 + 11) % 250);
    Require(stbi_write_png((folder / "reference.png").string().c_str(), 3, 2, 4, rgba.data(), 12) != 0,
            "PNG fixture write failed");
    auto negative_rgba = rgba; ++negative_rgba[1];
    Require(stbi_write_png((folder / "perturbed.png").string().c_str(), 3, 2, 4, negative_rgba.data(), 12) != 0,
            "PNG negative fixture write failed");
    const auto reference = Read(folder / "reference.png");
    Bytes metadata(reference.begin(), reference.begin() + 33);
    const std::string text("Forensic\0same pixels", 20);
    AppendBig32(metadata, uint32_t(text.size()));
    const size_t crc_start = metadata.size();
    metadata.insert(metadata.end(), {'t','E','X','t'});
    metadata.insert(metadata.end(), text.begin(), text.end());
    AppendBig32(metadata, Crc32(metadata.data() + crc_start, text.size() + 4));
    metadata.insert(metadata.end(), reference.begin() + 33, reference.end());
    Write(folder / "metadata.png", metadata);
    for (const std::string name : {"reference", "metadata", "perturbed"}) {
        auto payload = rgb;
        if (name == "perturbed") ++payload[1];
        const std::string header = name == "metadata" ? "P6\n# same pixels, extra metadata\n3 2\n255\n" : "P6\n3 2\n255\n";
        Bytes file(header.begin(), header.end()); file.insert(file.end(), payload.begin(), payload.end());
        Write(folder / (name + ".ppm"), file);
        ExrFixture(folder / (name + ".exr"), name == "metadata", name == "perturbed");
    }
    std::cout << "{\"fixtures_created\":true,\"width\":3,\"height\":2,\"formats\":[\"png\",\"exr\",\"ppm\"]}\n";
}

int Dimension(const std::string& value) {
    Require(!value.empty() && std::all_of(value.begin(), value.end(), [](unsigned char c) { return c >= '0' && c <= '9'; }), "invalid dimension argument");
    const auto parsed = std::stoull(value);
    Require(parsed > 0 && parsed <= 16384, "dimension argument outside domain");
    return int(parsed);
}

int main(int argc, char** argv) {
    try {
        if (argc == 3 && std::string(argv[1]) == "fixtures") { Fixtures(argv[2]); return 0; }
        Require(argc == 7 && std::string(argv[1]) == "compare", "usage: image_compare compare png|exr|ppm W H LEFT RIGHT; or fixtures NEW_DIRECTORY");
        return Compare(argv[2], Dimension(argv[3]), Dimension(argv[4]), argv[5], argv[6]);
    } catch (const std::exception& error) {
        std::cout << "{\"error\":" << Quoted(error.what()) << "}\n";
        return 2;
    }
}
