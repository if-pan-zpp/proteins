#pragma once
#include <google/protobuf/dynamic_message.h>
#include "line_parser.h"

namespace pdb::ir {
    class AutoLineParser: public LineParser {
    public:
        using Sectors = std::vector<std::pair<int, int>>;

    private:
        Sectors sectors;
        const google::protobuf::FieldDescriptor* header_in_anyfield;
        const google::protobuf::Descriptor* header_msg_type;
        google::protobuf::DynamicMessageFactory factory;

    public:
        AutoLineParser(Headers header, Sectors const& sectors);
        void parse(std::string const& line, AnyField *field) override;
    };
}
