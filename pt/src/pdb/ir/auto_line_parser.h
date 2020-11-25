#pragma once
#include "line_parser.h"

namespace pdb::ir {
    class AutoLineParser: public LineParser {
    public:
        using Sectors = std::vector<std::pair<int, int>>;

    private:
        Sectors sectors;
        const google::protobuf::FieldDescriptor* header_in_anyfield;
        const google::protobuf::Descriptor* header_msg_type;

    public:
        AutoLineParser(Headers header, Sectors const& sectors);
        void parse(std::string_view line, AnyField *field) override;
    };
}
