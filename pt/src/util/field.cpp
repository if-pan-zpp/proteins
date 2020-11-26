#include "field.h"
#include "strops.h"
using namespace std;
using namespace util;
using namespace google::protobuf;

Field::Field(Message *message,
        const FieldDescriptor *field,
        const Reflection *refl) {
    this->message = message;
    this->field = field;
    this->refl = refl ? refl : message->GetReflection();
}

void Field::setFrom(std::string_view s) {
    switch (field->cpp_type()) {
    case FieldDescriptor::CPPTYPE_INT32: {
        s = trim(s);
        if (s.empty()) break;

        auto value = stoi(string(s));
        refl->SetInt32(message, field, value);
        break;
    }
    case FieldDescriptor::CPPTYPE_STRING: {
        auto value = string(s);
        refl->SetString(message, field, move(value));
        break;
    }
    case FieldDescriptor::CPPTYPE_DOUBLE: {
        s = trim(s);
        if (s.empty()) break;

        auto value = stod(string(s));
        refl->SetDouble(message, field, value);
        break;
    }
    case FieldDescriptor::CPPTYPE_BOOL: {
        if (s == "T") {
            refl->SetBool(message, field, true);
        }
        else if (s == "F") {
            refl->SetBool(message, field, false);
        }
        break;
    }
    default: {}
    }
}
