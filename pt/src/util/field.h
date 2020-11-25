#pragma once
#include <google/protobuf/message.h>
#include <string>

namespace util {
    class Field {
    private:
        google::protobuf::Message* message;
        const google::protobuf::FieldDescriptor* field;
        const google::protobuf::Reflection* refl;

    public:
        Field(google::protobuf::Message* message,
            const google::protobuf::FieldDescriptor* fieldDesc,
            const google::protobuf::Reflection* reflection = nullptr);

        void setFrom(std::string_view s);
    };
}