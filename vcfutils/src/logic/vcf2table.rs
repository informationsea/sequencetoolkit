use crate::error::VCFUtilsError;
use crate::utils::tablewriter::{TableWriter, XlsxDataType, XlsxSheetWriter};
use nom::branch::alt;
use nom::bytes::complete::{tag, take_while1};
use nom::character::is_digit;
use nom::multi::many0;
use nom::sequence::tuple;
use std::collections::HashSet;
use std::fmt;
use std::io::{BufRead, Write};
use std::str;
use vcf::{self, U8Vec, VCFHeader, VCFReader, VCFRecord};

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct VCF2CSVConfig {
    pub split_multi_allelic: bool,
    pub decoded_genotype: bool,
    pub canonical_list: Option<HashSet<U8Vec>>,
    pub info_list: Vec<U8Vec>,
    pub format_list: Vec<U8Vec>,
    pub replace_sample_name: Option<Vec<U8Vec>>,
    pub group_names: Option<Vec<U8Vec>>,
}

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Clone, Copy)]
pub enum SnpEffImpact {
    High,
    Moderate,
    Low,
    Modifier,
}

impl SnpEffImpact {
    pub fn to_str(&self) -> &str {
        match self {
            SnpEffImpact::High => "HIGH",
            SnpEffImpact::Moderate => "MODERATE",
            SnpEffImpact::Low => "LOW",
            SnpEffImpact::Modifier => "MODIFIER",
        }
    }
}

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Clone)]
pub enum HeaderType {
    GroupName,
    VcfLine,
    AltIndex,
    GeneChange,
    GeneName,
    TranscriptName,
    AminoChange,
    CDSChange,
    CHROM,
    POS,
    ID,
    REF,
    ALT,
    QUAL,
    FILTER,
    SnpEff,
    SnpEffImpact(SnpEffImpact),
    Info(U8Vec, vcf::Number, vcf::ValueType, i32, String),
    // Genotype: SampleName, FormatID, Number, Type, Index, Description, ReplacedSampleName
    Genotype(
        U8Vec,
        U8Vec,
        vcf::Number,
        vcf::ValueType,
        i32,
        String,
        Option<U8Vec>,
    ),
    Empty,
}

fn add_number_suffix(
    s: &mut impl fmt::Write,
    num: &vcf::Number,
    index: i32,
) -> Result<(), fmt::Error> {
    match num {
        vcf::Number::Zero
        | vcf::Number::Number(1)
        | vcf::Number::Allele
        | vcf::Number::Other(_)
        | vcf::Number::Unknown
        | vcf::Number::Genotype => (),
        vcf::Number::Reference => match index {
            0 => {
                write!(s, "__Ref")?;
            }
            1 => {
                write!(s, "__Alt")?;
            }
            _ => unreachable!(),
        },
        vcf::Number::Number(_) => {
            write!(s, "__{}", index)?;
        }
    }
    Ok(())
}

impl ToString for HeaderType {
    fn to_string(&self) -> String {
        match self {
            HeaderType::GroupName => "Group Name".to_string(),
            HeaderType::VcfLine => "#".to_string(),
            HeaderType::AltIndex => "alt #".to_string(),
            HeaderType::GeneChange => "Gene change".to_string(),
            HeaderType::GeneName => "Gene".to_string(),
            HeaderType::TranscriptName => "Transcript".to_string(),
            HeaderType::AminoChange => "Amino change".to_string(),
            HeaderType::CDSChange => "CDS change".to_string(),
            HeaderType::CHROM => "CHROM".to_string(),
            HeaderType::POS => "POS".to_string(),
            HeaderType::ID => "ID".to_string(),
            HeaderType::REF => "REF".to_string(),
            HeaderType::ALT => "ALT".to_string(),
            HeaderType::QUAL => "QUAL".to_string(),
            HeaderType::FILTER => "FILTER".to_string(),
            HeaderType::Info(x, num, _, index, _) => {
                let mut s = str::from_utf8(x).unwrap().to_string();
                add_number_suffix(&mut s, num, *index).unwrap();
                s
            }
            HeaderType::Genotype(x, y, num, _, index, _, replace_sample_name) => {
                let mut s = if let Some(replace_sample_name) = replace_sample_name {
                    format!(
                        "{}__{}",
                        str::from_utf8(replace_sample_name).unwrap(),
                        str::from_utf8(y).unwrap()
                    )
                } else {
                    format!(
                        "{}__{}",
                        str::from_utf8(x).unwrap(),
                        str::from_utf8(y).unwrap()
                    )
                };
                add_number_suffix(&mut s, num, *index).unwrap();
                s
            }
            HeaderType::SnpEffImpact(impact) => format!("GeneImpact__{}", impact.to_str()),
            HeaderType::SnpEff => "SnpEff".to_string(),
            HeaderType::Empty => "".to_string(),
        }
    }
}

pub fn create_header_line(header: &VCFHeader, config: &VCF2CSVConfig) -> Vec<HeaderType> {
    let mut header_items = vec![
        HeaderType::VcfLine,
        HeaderType::CHROM,
        HeaderType::POS,
        HeaderType::ID,
        HeaderType::REF,
        HeaderType::ALT,
        HeaderType::QUAL,
        HeaderType::FILTER,
    ];

    if config.group_names.is_some() {
        header_items.insert(0, HeaderType::GroupName);
    }

    if config.split_multi_allelic {
        header_items.insert(1, HeaderType::AltIndex);
    }

    if header.info(b"ANN").is_some() {
        header_items.insert(2, HeaderType::SnpEff);
        header_items.insert(3, HeaderType::SnpEffImpact(SnpEffImpact::High));
        header_items.insert(4, HeaderType::SnpEffImpact(SnpEffImpact::Moderate));
        header_items.insert(5, HeaderType::SnpEffImpact(SnpEffImpact::Low));
        header_items.insert(6, HeaderType::SnpEffImpact(SnpEffImpact::Modifier));

        if config.canonical_list.is_some() {
            header_items.insert(2, HeaderType::GeneChange);
            header_items.insert(3, HeaderType::GeneName);
            header_items.insert(4, HeaderType::TranscriptName);
            header_items.insert(5, HeaderType::AminoChange);
            header_items.insert(6, HeaderType::CDSChange);
        }
    }

    for one_info in &config.info_list {
        let info = header.info(one_info).unwrap();
        match info.number {
            vcf::Number::Number(x) => {
                for i in 0..*x {
                    header_items.push(HeaderType::Info(
                        one_info.to_vec(),
                        info.number.clone(),
                        info.value_type.clone(),
                        i,
                        str::from_utf8(info.description).unwrap().to_string(),
                    ));
                }
            }
            vcf::Number::Reference => {
                header_items.push(HeaderType::Info(
                    one_info.to_vec(),
                    info.number.clone(),
                    info.value_type.clone(),
                    0,
                    str::from_utf8(info.description).unwrap().to_string(),
                ));
                header_items.push(HeaderType::Info(
                    one_info.to_vec(),
                    info.number.clone(),
                    info.value_type.clone(),
                    1,
                    str::from_utf8(info.description).unwrap().to_string(),
                ));
            }
            _ => {
                header_items.push(HeaderType::Info(
                    one_info.to_vec(),
                    info.number.clone(),
                    info.value_type.clone(),
                    0,
                    str::from_utf8(info.description).unwrap().to_string(),
                ));
            }
        }
    }

    for (sample_index, one_sample) in header.samples().iter().enumerate() {
        for one_format in &config.format_list {
            let format = header.format(one_format).unwrap();
            match format.number {
                vcf::Number::Number(x) => {
                    for i in 0..*x {
                        header_items.push(HeaderType::Genotype(
                            one_sample.to_vec(),
                            one_format.to_vec(),
                            format.number.clone(),
                            format.value_type.clone(),
                            i,
                            str::from_utf8(format.description).unwrap().to_string(),
                            config
                                .replace_sample_name
                                .as_ref()
                                .map(|x| x.get(sample_index).map(|y| y.to_vec()))
                                .flatten(),
                        ));
                    }
                }
                vcf::Number::Reference => {
                    header_items.push(HeaderType::Genotype(
                        one_sample.to_vec(),
                        one_format.to_vec(),
                        format.number.clone(),
                        format.value_type.clone(),
                        0,
                        str::from_utf8(format.description).unwrap().to_string(),
                        config
                            .replace_sample_name
                            .as_ref()
                            .map(|x| x.get(sample_index).map(|y| y.to_vec()))
                            .flatten(),
                    ));
                    header_items.push(HeaderType::Genotype(
                        one_sample.to_vec(),
                        one_format.to_vec(),
                        format.number.clone(),
                        format.value_type.clone(),
                        1,
                        str::from_utf8(format.description).unwrap().to_string(),
                        config
                            .replace_sample_name
                            .as_ref()
                            .map(|x| x.get(sample_index).map(|y| y.to_vec()))
                            .flatten(),
                    ));
                }
                _ => {
                    header_items.push(HeaderType::Genotype(
                        one_sample.to_vec(),
                        one_format.to_vec(),
                        format.number.clone(),
                        format.value_type.clone(),
                        0,
                        str::from_utf8(format.description).unwrap().to_string(),
                        config
                            .replace_sample_name
                            .as_ref()
                            .map(|x| x.get(sample_index).map(|y| y.to_vec()))
                            .flatten(),
                    ));
                }
            }
        }
    }

    header_items
}

fn write_comma_separated_values(writer: &mut U8Vec, values: &[U8Vec]) {
    for (i, d) in values.iter().enumerate() {
        if i != 0 {
            writer.push(b',');
        }
        writer.extend_from_slice(&d);
    }
}

fn write_value_for_alt_index(
    writer: &mut U8Vec,
    values: &[U8Vec],
    num: &vcf::Number,
    index: i32,
    alt_index: Option<usize>,
) {
    match num {
        vcf::Number::Zero
        | vcf::Number::Number(1)
        | vcf::Number::Other(_)
        | vcf::Number::Unknown
        | vcf::Number::Genotype => write_comma_separated_values(writer, values),
        vcf::Number::Allele => {
            if let Some(alt_index) = alt_index {
                if let Some(x) = values.get(alt_index) {
                    writer.extend_from_slice(x);
                }
            } else {
                write_comma_separated_values(writer, values);
            }
        }
        vcf::Number::Reference => match index {
            0 => {
                if let Some(x) = values.get(0) {
                    writer.extend_from_slice(x);
                }
            }
            1 => {
                if let Some(alt_index) = alt_index {
                    if let Some(x) = values.get(alt_index + 1) {
                        writer.extend_from_slice(x);
                    }
                } else {
                    write_comma_separated_values(writer, &values[1..])
                }
            }
            _ => unreachable!(),
        },
        vcf::Number::Number(_) => {
            if let Some(x) = values.get(index as usize) {
                writer.extend_from_slice(x);
            }
        }
    }
}

fn write_value_for_snpeff_ann(
    record: &VCFRecord,
    writer: &mut U8Vec,
    alt_index: usize,
) -> Result<(), VCFUtilsError> {
    let mut first_record = true;
    if let Some(snpeff_record) = record.info(b"ANN") {
        for one in snpeff_record {
            let mut ann = one.splitn(2, |x| *x == b'|');
            if ann.next() == Some(&record.alternative[alt_index]) {
                if first_record {
                    first_record = false;
                } else {
                    write!(writer, ",")?;
                }
                writer.write_all(one)?;
            }
        }
    }
    Ok(())
}

fn write_snpeff_all(
    record: &VCFRecord,
    writer: &mut U8Vec,
    alt_index: Option<usize>,
) -> Result<(), VCFUtilsError> {
    if let Some(snpeff_record) = record.info(b"ANN") {
        let alt_genotype: Option<&[u8]> = alt_index.map(|x| &record.alternative[x][..]);
        let mut first = true;
        for one in snpeff_record {
            let ann: Vec<_> = one.split(|x| *x == b'|').collect();
            let genotype = ann[0];
            if alt_genotype.is_none() || alt_genotype == Some(genotype) {
                if first {
                    first = false;
                } else {
                    writer.extend_from_slice(b", ");
                }
                writer.extend_from_slice(ann[6]);
                writer.extend_from_slice(b"(");
                writer.extend_from_slice(ann[3]);
                writer.extend_from_slice(b"):");
                if ann[10].is_empty() {
                    writer.extend_from_slice(ann[9]);
                } else {
                    writer.extend_from_slice(ann[10]);
                }
            }
        }
    }

    Ok(())
}

fn write_impact_gene(
    record: &VCFRecord,
    writer: &mut U8Vec,
    alt_index: Option<usize>,
    header: &HeaderType,
) -> Result<(), VCFUtilsError> {
    if let Some(snpeff_record) = record.info(b"ANN") {
        let alt_genotype: Option<&[u8]> = alt_index.map(|x| &record.alternative[x][..]);
        let mut genes: HashSet<U8Vec> = HashSet::new();
        for one in snpeff_record {
            let ann: Vec<_> = one.split(|x| *x == b'|').collect();
            let genotype = ann[0];
            if alt_genotype.is_none() || alt_genotype == Some(genotype) {
                match header {
                    HeaderType::SnpEffImpact(impact) => {
                        if ann[2] == impact.to_str().as_bytes() {
                            genes.insert(ann[3].to_vec());
                        }
                    }
                    _ => unreachable!(),
                }
            }
        }
        let mut gene_list: Vec<_> = genes.iter().collect();
        gene_list.sort();
        for (i, one) in gene_list.iter().enumerate() {
            if i != 0 {
                writer.extend_from_slice(b", ");
            }
            writer.extend_from_slice(one);
        }
    }
    Ok(())
}

fn write_canonical_gene(
    record: &VCFRecord,
    writer: &mut U8Vec,
    alt_index: Option<usize>,
    header: &HeaderType,
    canonical_list: Option<&HashSet<U8Vec>>,
) -> Result<(), VCFUtilsError> {
    if let Some(snpeff_record) = record.info(b"ANN") {
        let alt_genotype: Option<&[u8]> = alt_index.map(|x| &record.alternative[x][..]);
        let mut first = true;
        for one in snpeff_record {
            let ann: Vec<_> = one.split(|x| *x == b'|').collect();
            let genotype = ann[0];
            let transcript = ann[6];
            if (alt_genotype.is_none() || alt_genotype == Some(genotype))
                && canonical_list
                    .map(|x| x.contains(transcript))
                    .unwrap_or(true)
            {
                if first {
                    first = false;
                } else {
                    write!(writer, ", ")?;
                }
                match header {
                    HeaderType::GeneChange => {
                        write!(writer, "{} ", str::from_utf8(ann[3])?)?;
                        if ann[10].is_empty() {
                            write!(writer, "{}", str::from_utf8(ann[9])?)?;
                        } else {
                            write!(writer, "{}", str::from_utf8(ann[10])?)?
                        }
                    }
                    HeaderType::GeneName => write!(writer, "{}", str::from_utf8(ann[3])?)?,
                    HeaderType::TranscriptName => write!(writer, "{}", str::from_utf8(ann[6])?)?,
                    HeaderType::AminoChange => write!(writer, "{}", str::from_utf8(ann[10])?)?,
                    HeaderType::CDSChange => write!(writer, "{}", str::from_utf8(ann[9])?)?,
                    _ => unreachable!(),
                }
            }
        }
    }
    Ok(())
}

fn genotype_for_index(record: &VCFRecord, index: usize) -> &[u8] {
    if index == 0 {
        &record.reference[..]
    } else {
        record
            .alternative
            .get(index - 1)
            .map(|x| &x[..])
            .unwrap_or(&b"?"[..])
    }
}

fn write_value_for_decoded_genotype(
    writer: &mut U8Vec,
    record: &VCFRecord,
    genotype_record: &[U8Vec],
) -> Result<(), VCFUtilsError> {
    if let Some(one) = genotype_record.get(0) {
        if let Ok((_, parsed)) = tuple::<_, _, (&[u8], _), _>((
            many0(tuple((
                take_while1(|x| is_digit(x) || x == b'.'),
                alt((tag(b"|"), tag(b"/"))),
            ))),
            take_while1(|x| is_digit(x) || x == b'.'),
        ))(&one[..])
        {
            for x in parsed.0 {
                if x.0 == b"." {
                    writer.write_all(b".")?;
                } else {
                    writer.write_all(genotype_for_index(
                        record,
                        str::from_utf8(x.0)?.parse::<usize>()?,
                    ))?;
                }
                writer.write_all(x.1)?;
            }
            if parsed.1 == b"." {
                writer.write_all(b".")?;
            } else {
                writer.write_all(genotype_for_index(
                    record,
                    str::from_utf8(parsed.1)?.parse::<usize>()?,
                ))?;
            }
        }
    }
    Ok(())
}

pub fn vcf2table_set_data_type(
    header_contents: &[HeaderType],
    writer: &mut XlsxSheetWriter,
) -> Result<(), VCFUtilsError> {
    let types: Vec<_> = header_contents
        .iter()
        .map(|x| match x {
            HeaderType::POS | HeaderType::QUAL | HeaderType::VcfLine | HeaderType::AltIndex => {
                XlsxDataType::Number
            }
            HeaderType::Info(_, _, t, _, _) => match t {
                vcf::ValueType::Float | vcf::ValueType::Integer => XlsxDataType::Number,
                _ => XlsxDataType::String,
            },
            HeaderType::Genotype(_, _, _, t, _, _, _) => match t {
                vcf::ValueType::Float | vcf::ValueType::Integer => XlsxDataType::Number,
                _ => XlsxDataType::String,
            },
            _ => XlsxDataType::String,
        })
        .collect();
    let comments: Vec<_> = header_contents
        .iter()
        .map(|x| match x {
            HeaderType::POS | HeaderType::QUAL | HeaderType::VcfLine => "".to_string(),
            HeaderType::Info(_, _, _, _, comment) => comment.to_string(),
            HeaderType::Genotype(_, _, _, _, _, comment, _) => comment.to_string(),
            HeaderType::SnpEffImpact(impact) => {
                format!("{} impact (See \"Effect prediction details\" at http://snpeff.sourceforge.net/SnpEff_manual.html)", impact.to_str())
            }
            HeaderType::SnpEff => {
                "All transcript annotation by snpEff".to_string()
            }
            _ => "".to_string(),
        })
        .collect();
    writer.set_data_type(&types);
    writer.set_header_comment(&comments);
    Ok(())
}

fn setup_row(
    group_name: Option<&U8Vec>,
    header_contents: &[HeaderType],
    record: &VCFRecord,
    row: &mut Vec<U8Vec>,
    index: u32,
    alt_index: Option<usize>,
    translate_genotype: bool,
    canonical_list: Option<&HashSet<U8Vec>>,
) -> Result<(), VCFUtilsError> {
    for (header, column) in header_contents.iter().zip(row.iter_mut()) {
        column.clear();
        match header {
            HeaderType::GroupName => {
                if let Some(g) = group_name {
                    column.extend_from_slice(g);
                }
            }
            HeaderType::VcfLine => {
                write!(column, "{}", index)?;
            }
            HeaderType::AltIndex => {
                if let Some(alt_index) = alt_index {
                    write!(column, "{}", alt_index + 1)?;
                }
            }
            HeaderType::GeneChange => {
                write_canonical_gene(record, column, alt_index, header, canonical_list)?
            }
            HeaderType::GeneName => {
                write_canonical_gene(record, column, alt_index, header, canonical_list)?
            }
            HeaderType::TranscriptName => {
                write_canonical_gene(record, column, alt_index, header, canonical_list)?
            }
            HeaderType::AminoChange => {
                write_canonical_gene(record, column, alt_index, header, canonical_list)?
            }
            HeaderType::CDSChange => {
                write_canonical_gene(record, column, alt_index, header, canonical_list)?
            }
            HeaderType::CHROM => column.extend_from_slice(&record.chromosome),
            HeaderType::POS => {
                write!(column, "{}", record.position)?;
            }
            HeaderType::ID => {
                write_comma_separated_values(column, &record.id);
            }
            HeaderType::REF => column.extend_from_slice(&record.reference),
            HeaderType::ALT => {
                if let Some(alt_index) = alt_index {
                    column.extend_from_slice(&record.alternative[alt_index]);
                } else {
                    write_comma_separated_values(column, &record.alternative);
                }
            }
            HeaderType::QUAL => {
                if let Some(q) = record.qual {
                    write!(column, "{}", q)?;
                }
            }
            HeaderType::FILTER => {
                write_comma_separated_values(column, &record.filter);
            }
            HeaderType::Info(key, number, value_type, index, _) => {
                if value_type == &vcf::ValueType::Flag {
                    if record.info(key).is_some() {
                        write!(column, "TRUE")?;
                    } else {
                        write!(column, "FALSE")?;
                    }
                } else if let Some(values) = record.info(key) {
                    if key == b"ANN" && number == &vcf::Number::Unknown {
                        if let Some(alt_index) = alt_index {
                            write_value_for_snpeff_ann(record, column, alt_index)?;
                        } else {
                            write_value_for_alt_index(column, values, number, *index, alt_index);
                        }
                    } else {
                        write_value_for_alt_index(column, values, number, *index, alt_index);
                    }
                }
            }
            HeaderType::Genotype(sample_name, key, number, _, index, _, _) => {
                if let Some(values) = record.genotype(sample_name, key) {
                    if translate_genotype && (key == b"GT" || key == b"PGT") {
                        write_value_for_decoded_genotype(column, record, values)?;
                    } else {
                        write_value_for_alt_index(column, values, number, *index, alt_index);
                    }
                }
            }
            HeaderType::SnpEffImpact(_) => {
                write_impact_gene(record, column, alt_index, header)?;
            }
            HeaderType::SnpEff => {
                write_snpeff_all(record, column, alt_index)?;
            }
            HeaderType::Empty => {}
        }
    }

    Ok(())
}

pub fn vcf2table<R: BufRead, W: TableWriter>(
    vcf_reader: &mut VCFReader<R>,
    header_contents: &[HeaderType],
    config: &VCF2CSVConfig,
    group_name: Option<&U8Vec>,
    write_header: bool,
    mut writer: W,
) -> Result<u32, VCFUtilsError> {
    let header_string: Vec<_> = header_contents.iter().map(|x| x.to_string()).collect();
    writer.set_header(&header_string);
    if write_header {
        writer.write_header()?;
    }
    let mut row: Vec<U8Vec> = header_contents.iter().map(|_| Vec::new()).collect();

    let mut index: u32 = 0;

    let mut record = VCFRecord::new(vcf_reader.header().clone());
    while vcf_reader.next_record(&mut record)? {
        index += 1;

        if config.split_multi_allelic {
            for (alt_index, _) in record.alternative.iter().enumerate() {
                setup_row(
                    group_name,
                    header_contents,
                    &record,
                    &mut row,
                    index,
                    Some(alt_index),
                    config.decoded_genotype,
                    config.canonical_list.as_ref(),
                )?;
                writer.write_row_bytes(&row.iter().map(|x| -> &[u8] { &x }).collect::<Vec<_>>())?;
            }
        } else {
            setup_row(
                group_name,
                header_contents,
                &record,
                &mut row,
                index,
                None,
                config.decoded_genotype,
                config.canonical_list.as_ref(),
            )?;
            writer.write_row_bytes(&row.iter().map(|x| -> &[u8] { &x }).collect::<Vec<_>>())?;
        }
    }

    Ok(index)
}

pub fn merge_header_contents(original: &[HeaderType], new: &[HeaderType]) -> Vec<HeaderType> {
    let mut merged = Vec::new();

    for one in original {
        if let Some(found) = new
            .iter()
            .filter(|x| match *x {
                HeaderType::Genotype(
                    x_sample_name,
                    x_format_id,
                    _,
                    _,
                    x_index,
                    _,
                    x_replace_sample_name,
                ) => {
                    if let HeaderType::Genotype(
                        one_sample_name,
                        one_format_id,
                        _,
                        _,
                        one_index,
                        _,
                        one_replace_sample_name,
                    ) = one
                    {
                        if x_format_id != one_format_id || x_index != one_index {
                            return false;
                        }
                        if let Some(x_replace_sample_name) = x_replace_sample_name {
                            if let Some(one_replace_sample_name) = one_replace_sample_name {
                                return x_replace_sample_name == one_replace_sample_name;
                            }
                        }
                        one_sample_name == x_sample_name
                    } else {
                        false
                    }
                }
                HeaderType::Info(x_id, _, _, x_index, _) => {
                    if let HeaderType::Info(one_id, _, _, one_index, _) = one {
                        one_id == x_id && one_index == x_index
                    } else {
                        false
                    }
                }
                y => y == one,
            })
            .next()
        {
            merged.push(found.clone());
        } else {
            merged.push(HeaderType::Empty);
        }
    }

    merged
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::utils::tablewriter;
    use std::io::BufReader;

    #[test]
    fn test_write_comma_separated_values() {
        let mut result = Vec::new();
        write_comma_separated_values(&mut result, &[b"AB".to_vec(), b"CD".to_vec()][..]);
        assert_eq!(result, b"AB,CD");
        result.clear();
        write_comma_separated_values(&mut result, &[b"AB".to_vec()][..]);
        assert_eq!(result, b"AB");
        result.clear();
        write_comma_separated_values(
            &mut result,
            &[b"AB".to_vec(), b"CD".to_vec(), b"EF".to_vec()][..],
        );
        assert_eq!(result, b"AB,CD,EF");
    }

    #[test]
    fn test_write_value_for_alt_index() {
        let mut buffer = Vec::new();
        let values = &[b"A".to_vec(), b"B".to_vec(), b"C".to_vec()][..];
        write_value_for_alt_index(&mut buffer, values, &vcf::Number::Unknown, 0, Some(0));
        assert_eq!(buffer, b"A,B,C");
        buffer.clear();

        write_value_for_alt_index(&mut buffer, values, &vcf::Number::Number(3), 0, Some(0));
        assert_eq!(buffer, b"A");
        buffer.clear();

        write_value_for_alt_index(&mut buffer, values, &vcf::Number::Number(3), 1, Some(0));
        assert_eq!(buffer, b"B");
        buffer.clear();

        write_value_for_alt_index(&mut buffer, values, &vcf::Number::Number(3), 2, Some(0));
        assert_eq!(buffer, b"C");
        buffer.clear();

        write_value_for_alt_index(&mut buffer, values, &vcf::Number::Allele, 0, Some(0));
        assert_eq!(buffer, b"A");
        buffer.clear();

        write_value_for_alt_index(&mut buffer, values, &vcf::Number::Allele, 0, Some(1));
        assert_eq!(buffer, b"B");
        buffer.clear();

        write_value_for_alt_index(&mut buffer, values, &vcf::Number::Reference, 0, Some(1));
        assert_eq!(buffer, b"A");
        buffer.clear();

        write_value_for_alt_index(&mut buffer, values, &vcf::Number::Reference, 1, Some(1));
        assert_eq!(buffer, b"C");
        buffer.clear();
    }

    #[test]
    fn test_vcf2table_csv_split_multi() -> Result<(), VCFUtilsError> {
        let vcf_data = include_bytes!("../../testfiles/simple1.vcf");
        let config = VCF2CSVConfig {
            split_multi_allelic: true,
            decoded_genotype: false,
            canonical_list: None,
            info_list: vec![
                b"AC".to_vec(),
                b"AF".to_vec(),
                b"AN".to_vec(),
                b"DP".to_vec(),
            ],
            format_list: vec![b"AD".to_vec(), b"DP".to_vec(), b"GT".to_vec()],
            replace_sample_name: None,
            group_names: None,
        };
        let mut vcf_data_reader = BufReader::new(&vcf_data[..]);
        let mut vcf_reader = vcf::VCFReader::new(&mut vcf_data_reader)?;
        let mut write_bytes = Vec::<u8>::new();
        let header_contents = create_header_line(&vcf_reader.header(), &config);
        vcf2table(
            &mut vcf_reader,
            &header_contents,
            &config,
            None,
            true,
            &mut tablewriter::CSVWriter::new(&mut write_bytes),
        )?;
        std::fs::File::create("../target/simple1-expected-multiallelic-split.csv")?
            .write_all(&write_bytes)?;
        assert_eq!(
            &write_bytes[..],
            &include_bytes!("../../testfiles/simple1-expected-multiallelic-split.csv")[..]
        );
        Ok(())
    }

    #[test]
    fn test_vcf2table_csv_split_multi_snpeff() -> Result<(), VCFUtilsError> {
        let vcf_data = include_bytes!("../../testfiles/simple1-snpeff.vcf");
        let config = VCF2CSVConfig {
            split_multi_allelic: true,
            decoded_genotype: true,
            canonical_list: None,
            info_list: vec![
                b"AC".to_vec(),
                b"AF".to_vec(),
                b"AN".to_vec(),
                b"DP".to_vec(),
                b"ANN".to_vec(),
                b"FLAG".to_vec(),
            ],
            format_list: vec![b"AD".to_vec(), b"DP".to_vec(), b"GT".to_vec()],
            replace_sample_name: None,
            group_names: None,
        };
        let mut vcf_data_reader = BufReader::new(&vcf_data[..]);
        let mut vcf_reader = vcf::VCFReader::new(&mut vcf_data_reader)?;
        let mut write_bytes = Vec::<u8>::new();
        let header_contents = create_header_line(&vcf_reader.header(), &config);
        vcf2table(
            &mut vcf_reader,
            &header_contents,
            &config,
            None,
            true,
            &mut tablewriter::CSVWriter::new(&mut write_bytes),
        )?;
        std::fs::File::create("../target/split-snpeff.csv")?.write_all(&write_bytes)?;
        assert_eq!(
            &write_bytes[..],
            &include_bytes!("../../testfiles/simple1-expected-multiallelic-split-snpeff.csv")[..]
        );
        Ok(())
    }

    #[test]
    fn test_vcf2table_csv_split_multi_snpeff_with_canonical() -> Result<(), VCFUtilsError> {
        let vcf_data = include_bytes!("../../testfiles/1kGP-subset-snpeff.vcf");
        let config = VCF2CSVConfig {
            split_multi_allelic: true,
            decoded_genotype: true,
            canonical_list: Some(
                [b"ENST00000380152.7_1"]
                    .iter()
                    .map(|x| x.to_vec())
                    .collect(),
            ),
            info_list: vec![
                b"AC".to_vec(),
                b"AF".to_vec(),
                b"AN".to_vec(),
                b"DP".to_vec(),
                b"ANN".to_vec(),
            ],
            format_list: vec![b"AD".to_vec(), b"DP".to_vec(), b"GT".to_vec()],
            replace_sample_name: None,
            group_names: None,
        };
        let mut vcf_data_reader = BufReader::new(&vcf_data[..]);
        let mut vcf_reader = vcf::VCFReader::new(&mut vcf_data_reader)?;
        let mut write_bytes = Vec::<u8>::new();
        let header_contents = create_header_line(&vcf_reader.header(), &config);
        vcf2table(
            &mut vcf_reader,
            &header_contents,
            &config,
            None,
            true,
            &mut tablewriter::CSVWriter::new(&mut write_bytes),
        )?;
        std::fs::File::create("../target/split-snpeff-canonical.csv")?.write_all(&write_bytes)?;
        // assert_eq!(
        //     &write_bytes[..],
        //     &include_bytes!("../../testfiles/simple1-expected-multiallelic-split-snpeff.csv")[..]
        // );
        Ok(())
    }

    #[test]
    fn test_vcf2table_csv_no_split() -> Result<(), VCFUtilsError> {
        let vcf_data = include_bytes!("../../testfiles/simple1.vcf");
        let config = VCF2CSVConfig {
            split_multi_allelic: false,
            decoded_genotype: false,
            canonical_list: None,
            info_list: vec![
                b"AC".to_vec(),
                b"AF".to_vec(),
                b"AN".to_vec(),
                b"DP".to_vec(),
            ],
            format_list: vec![b"AD".to_vec(), b"DP".to_vec(), b"GT".to_vec()],
            replace_sample_name: None,
            group_names: None,
        };
        let mut vcf_data_reader = BufReader::new(&vcf_data[..]);
        let mut vcf_reader = vcf::VCFReader::new(&mut vcf_data_reader)?;
        let mut write_bytes = Vec::<u8>::new();
        let header_contents = create_header_line(&vcf_reader.header(), &config);
        vcf2table(
            &mut vcf_reader,
            &header_contents,
            &config,
            None,
            true,
            &mut tablewriter::CSVWriter::new(&mut write_bytes),
        )?;
        assert_eq!(
            &write_bytes[..],
            &include_bytes!("../../testfiles/simple1-expected-no-split.csv")[..]
        );
        Ok(())
    }

    #[test]
    fn test_vcf2table_xlsx_split_multi() -> Result<(), VCFUtilsError> {
        let vcf_data = include_bytes!("../../testfiles/simple1.vcf");
        let config = VCF2CSVConfig {
            split_multi_allelic: true,
            decoded_genotype: false,
            canonical_list: None,
            info_list: vec![
                b"AC".to_vec(),
                b"AN".to_vec(),
                b"AF".to_vec(),
                b"DP".to_vec(),
            ],
            format_list: vec![b"GT".to_vec(), b"AD".to_vec(), b"DP".to_vec()],
            replace_sample_name: None,
            group_names: None,
        };
        let mut vcf_data_reader = BufReader::new(&vcf_data[..]);
        let mut vcf_reader = vcf::VCFReader::new(&mut vcf_data_reader)?;
        let workbook = xlsxwriter::Workbook::new("../target/table-split-multi.xlsx");
        let mut sheet = workbook.add_worksheet(None)?;
        let mut writer = tablewriter::XlsxSheetWriter::new(&mut sheet);
        let header_contents = create_header_line(&vcf_reader.header(), &config);
        vcf2table_set_data_type(&header_contents, &mut writer)?;
        let row = vcf2table(
            &mut vcf_reader,
            &header_contents,
            &config,
            None,
            true,
            &mut writer,
        )?;
        sheet.autofilter(0, 0, row, (header_contents.len() - 1) as u16)?;
        workbook.close()?;
        Ok(())
    }

    #[test]
    fn test_vcf2table_xlsx_no_split_multi() -> Result<(), VCFUtilsError> {
        let vcf_data = include_bytes!("../../testfiles/simple1.vcf");
        let config = VCF2CSVConfig {
            split_multi_allelic: false,
            decoded_genotype: false,
            canonical_list: None,
            info_list: vec![
                b"AC".to_vec(),
                b"AN".to_vec(),
                b"AF".to_vec(),
                b"DP".to_vec(),
            ],
            format_list: vec![b"GT".to_vec(), b"AD".to_vec(), b"DP".to_vec()],
            replace_sample_name: None,
            group_names: None,
        };
        let mut vcf_data_reader = BufReader::new(&vcf_data[..]);
        let mut vcf_reader = vcf::VCFReader::new(&mut vcf_data_reader)?;
        let workbook = xlsxwriter::Workbook::new("../target/table-no-split.xlsx");
        let mut sheet = workbook.add_worksheet(None)?;
        let mut writer = tablewriter::XlsxSheetWriter::new(&mut sheet);
        let header_contents = create_header_line(&vcf_reader.header(), &config);
        vcf2table_set_data_type(&header_contents, &mut writer)?;
        let row = vcf2table(
            &mut vcf_reader,
            &header_contents,
            &config,
            None,
            true,
            &mut writer,
        )?;
        sheet.autofilter(0, 0, row, (header_contents.len() - 1) as u16)?;
        workbook.close()?;
        Ok(())
    }

    #[test]
    fn test_vcf2table_xlsx_split_multi_with_group_name() -> Result<(), VCFUtilsError> {
        let vcf_data = include_bytes!("../../testfiles/simple1.vcf");
        let config = VCF2CSVConfig {
            split_multi_allelic: true,
            decoded_genotype: false,
            canonical_list: None,
            info_list: vec![
                b"AC".to_vec(),
                b"AN".to_vec(),
                b"AF".to_vec(),
                b"DP".to_vec(),
            ],
            format_list: vec![b"GT".to_vec(), b"AD".to_vec(), b"DP".to_vec()],
            replace_sample_name: Some(vec![b"SAMPLE1".to_vec()]),
            group_names: Some(vec![b"GROUP".to_vec()]),
        };
        let mut vcf_data_reader = BufReader::new(&vcf_data[..]);
        let mut vcf_reader = vcf::VCFReader::new(&mut vcf_data_reader)?;
        let workbook = xlsxwriter::Workbook::new("../target/table-split-multi-with-group.xlsx");
        let mut sheet = workbook.add_worksheet(None)?;
        let mut writer = tablewriter::XlsxSheetWriter::new(&mut sheet);
        let header_contents = create_header_line(&vcf_reader.header(), &config);
        vcf2table_set_data_type(&header_contents, &mut writer)?;
        let row = vcf2table(
            &mut vcf_reader,
            &header_contents,
            &config,
            None,
            true,
            &mut writer,
        )?;
        sheet.autofilter(0, 0, row, (header_contents.len() - 1) as u16)?;
        workbook.close()?;
        Ok(())
    }
}
