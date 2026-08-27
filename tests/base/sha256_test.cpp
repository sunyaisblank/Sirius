#include "sirius/base/sha256.h"

#include <gtest/gtest.h>

#include "support/scoped_temporary_directory.h"

#include <fstream>
#include <string>
#include <string_view>

namespace sirius::base {
namespace {

TEST(Sha256Tests, EmptyStringMatchesNistKnownAnswer) {
    const auto digest = Sha256Hex(std::string_view{});
    ASSERT_TRUE(digest.has_value()) << digest.error();
    EXPECT_EQ(*digest,
              "e3b0c44298fc1c149afbf4c8996fb924"
              "27ae41e4649b934ca495991b7852b855");
}

TEST(Sha256Tests, AbcMatchesNistKnownAnswer) {
    const auto digest = Sha256Hex("abc");
    ASSERT_TRUE(digest.has_value()) << digest.error();
    EXPECT_EQ(*digest,
              "ba7816bf8f01cfea414140de5dae2223"
              "b00361a396177a9cb410ff61f20015ad");
}

TEST(Sha256Tests, MultiBlockPaddingMatchesNistKnownAnswer) {
    const auto digest = Sha256Hex(
        "abcdbcdecdefdefgefghfghighijhijk"
        "ijkljklmklmnlmnomnopnopq");
    ASSERT_TRUE(digest.has_value()) << digest.error();
    EXPECT_EQ(*digest,
              "248d6a61d20638b8e5c026930c3e6039"
              "a33ce45964ff2167f6ecedd419db06c1");
}

TEST(Sha256Tests, FileStreamingMatchesTheSameKnownAnswer) {
    sirius::test::ScopedTemporaryDirectory temporary("sirius-sha256");
    const auto path = temporary.path() / "abc.bin";
    {
        std::ofstream output(path, std::ios::binary);
        ASSERT_TRUE(output.good());
        output << "abc";
    }
    const auto digest = Sha256File(path);
    ASSERT_TRUE(digest.has_value()) << digest.error();
    EXPECT_EQ(*digest,
              "ba7816bf8f01cfea414140de5dae2223"
              "b00361a396177a9cb410ff61f20015ad");
}

TEST(Sha256Tests, MillionByteFileMatchesNistKnownAnswerAcrossReadChunks) {
    const std::string payload(1'000'000U, 'a');
    const auto memory_digest = Sha256Hex(payload);
    ASSERT_TRUE(memory_digest.has_value()) << memory_digest.error();
    EXPECT_EQ(*memory_digest,
              "cdc76e5c9914fb9281a1c7e284d73e67"
              "f1809a48a497200e046d39ccc7112cd0");

    sirius::test::ScopedTemporaryDirectory temporary("sirius-sha256-million");
    const auto path = temporary.path() / "million-a.bin";
    {
        std::ofstream output(path, std::ios::binary);
        ASSERT_TRUE(output.good());
        output.write(payload.data(), static_cast<std::streamsize>(payload.size()));
        ASSERT_TRUE(output.good());
    }
    const auto file_digest = Sha256File(path);
    ASSERT_TRUE(file_digest.has_value()) << file_digest.error();
    EXPECT_EQ(*file_digest, *memory_digest);
}

}  // namespace
}  // namespace sirius::base
