#include "auto_line_parser.h"
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

void pdb::ir::AutoLineParser::parse(const std::string &line,
        pdb::ir::AnyField *field) {
    auto field_value_prot = factory.GetPrototype(header_msg_type);
    auto field_value = field_value_prot->New();
    auto fv_refl = field_value->GetReflection();

    field->GetReflection()->SetAllocatedMessage(
        field, field_value, header_in_anyfield);
    for (int i = 0; i < header_msg_type->field_count(); ++i) {
        auto header_field_type = header_msg_type->field(i);
        auto [start, end] = sectors[i];
        auto sector = line.substr(start, end - start);

        switch (header_field_type->cpp_type()) {
        case FieldDescriptor::CPPTYPE_INT32: {
            int32_t value = stoi(sector);
            fv_refl->SetInt32(field_value, header_field_type, value);
            break;
        }
        case FieldDescriptor::CPPTYPE_DOUBLE: {
            double value = stod(sector);
            fv_refl->SetDouble(field_value, header_field_type, value);
            break;
        }
        case FieldDescriptor::CPPTYPE_STRING: {
            fv_refl->SetString(field_value, header_field_type, sector);
            break;
        }
        default: {}
        }
    }
}
