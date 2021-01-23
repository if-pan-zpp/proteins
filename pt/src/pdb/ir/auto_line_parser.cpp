#include "auto_line_parser.hpp"
#include "util/field.hpp"
using namespace google::protobuf;
using namespace std;

pdb::ir::AutoLineParser::AutoLineParser(pdb::ir::Headers header,
        const pdb::ir::AutoLineParser::Sectors &sectors) {
    this->sectors = sectors;
    auto field_desc = Field::GetDescriptor();
    header_in_anyfield = field_desc->FindFieldByNumber(1 + header);
    header_msg_type = header_in_anyfield->message_type();

    assert(sectors.size() == header_msg_type->field_count());
}

void pdb::ir::AutoLineParser::parse(std::string_view line,
        pdb::ir::Field *field) {
    auto field_refl = Field::GetReflection();
    auto field_value = field_refl->MutableMessage(field, header_in_anyfield);
    auto field_value_refl = field_value->GetReflection();

    for (int i = 0; i < header_msg_type->field_count(); ++i) {
        auto header_field_type = header_msg_type->field(i);
        auto [start, end] = sectors[i];
        string_view sector(line.data() + start - 1, end - start + 1);

        auto ufield = ::util::Field(field_value, header_field_type, field_value_refl);
        ufield.setFrom(sector);
    }
}
