#include "auto_line_parser.h"
#include "util/field.h"
using namespace google::protobuf;
using namespace std;

pdb::ir::AutoLineParser::AutoLineParser(pdb::ir::Headers header,
        const pdb::ir::AutoLineParser::Sectors &sectors) {
    this->sectors = sectors;
    auto af_desc = AnyField::GetDescriptor();
    header_in_anyfield = af_desc->FindFieldByNumber(1 + header);
    header_msg_type = header_in_anyfield->message_type();

    assert(sectors.size() == header_msg_type->field_count());
}

void pdb::ir::AutoLineParser::parse(std::string_view line,
        pdb::ir::AnyField *field) {
    auto field_refl = AnyField::GetReflection();
    auto af_value = field_refl->MutableMessage(field, header_in_anyfield);
    auto af_value_refl = af_value->GetReflection();

    for (int i = 0; i < header_msg_type->field_count(); ++i) {
        auto header_field_type = header_msg_type->field(i);
        auto [start, end] = sectors[i];
        string_view sector(line.data() + start - 1, end - start + 1);

        auto ufield = ::util::Field(af_value, header_field_type, af_value_refl);
        ufield.setFrom(sector);
    }
}
