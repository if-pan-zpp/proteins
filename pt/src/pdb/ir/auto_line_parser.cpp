#include "auto_line_parser.h"
#include "util/trim.h"
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
    auto field_refl = field->GetReflection();
    auto af_value = field_refl->MutableMessage(field, header_in_anyfield);
    auto af_value_refl = af_value->GetReflection();

    for (int i = 0; i < header_msg_type->field_count(); ++i) {
        auto header_field_type = header_msg_type->field(i);
        auto [start, end] = sectors[i];
        auto sector = line.substr(start, end - start);

        switch (header_field_type->cpp_type()) {
        case FieldDescriptor::CPPTYPE_INT32: {
            sector = ::util::trim(sector);
            if (sector.empty()) break;

            int32_t value = stoi(sector);
            af_value_refl->SetInt32(af_value, header_field_type, value);
            break;
        }
        case FieldDescriptor::CPPTYPE_DOUBLE: {
            sector = ::util::trim(sector);
            if (sector.empty()) break;

            double value = stod(sector);
            af_value_refl->SetDouble(af_value, header_field_type, value);
            break;
        }
        case FieldDescriptor::CPPTYPE_STRING: {
            af_value_refl->SetString(af_value, header_field_type, sector);
            break;
        }
        default: {}
        }
    }
}
