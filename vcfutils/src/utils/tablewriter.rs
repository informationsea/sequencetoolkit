use crate::error::VCFUtilsError;
use std::collections::HashMap;
use std::io::Write;
use std::str;

pub trait TableWriter {
    fn set_header(&mut self, items: &[String]);
    fn header(&self) -> &[String];
    fn write_row(&mut self, items: &[&str]) -> Result<(), VCFUtilsError>;
    fn write_header(&mut self) -> Result<(), VCFUtilsError> {
        let header: Vec<_> = self.header().iter().map(|x| x.to_string()).collect();
        self.write_row(&header.iter().map(|x| x.as_str()).collect::<Vec<&str>>())
    }

    fn write_row_bytes(&mut self, items: &[&[u8]]) -> Result<(), VCFUtilsError> {
        self.write_row(&items.iter().try_fold::<_, _, Result<_, VCFUtilsError>>(
            Vec::new(),
            |mut v, x| {
                let s = str::from_utf8(x)?;
                v.push(s);
                Ok(v)
            },
        )?)
    }

    fn write_dict(&mut self, items: &HashMap<&str, &str>) -> Result<(), VCFUtilsError> {
        let mut row: Vec<&str> = Vec::new();
        for one in self.header() {
            let s: &str = one;
            row.push(items.get(s).copied().unwrap_or(""));
        }
        self.write_row(&row)
    }

    fn write_dict_bytes(&mut self, items: &HashMap<&[u8], &[u8]>) -> Result<(), VCFUtilsError> {
        let mut row: Vec<&[u8]> = Vec::new();
        for one in self.header() {
            row.push(items.get(one.as_bytes()).copied().unwrap_or(b""));
        }
        self.write_row_bytes(&row)
    }
}

impl<T: TableWriter + ?Sized> TableWriter for Box<T> {
    fn set_header(&mut self, items: &[String]) {
        (**self).set_header(items)
    }
    fn header(&self) -> &[String] {
        (**self).header()
    }
    fn write_row(&mut self, items: &[&str]) -> Result<(), VCFUtilsError> {
        (**self).write_row(items)
    }
    fn write_header(&mut self) -> Result<(), VCFUtilsError> {
        (**self).write_header()
    }
    fn write_row_bytes(&mut self, items: &[&[u8]]) -> Result<(), VCFUtilsError> {
        (**self).write_row_bytes(items)
    }
    fn write_dict(&mut self, items: &HashMap<&str, &str>) -> Result<(), VCFUtilsError> {
        (**self).write_dict(items)
    }
    fn write_dict_bytes(&mut self, items: &HashMap<&[u8], &[u8]>) -> Result<(), VCFUtilsError> {
        (**self).write_dict_bytes(items)
    }
}

impl<T: TableWriter + ?Sized> TableWriter for &mut T {
    fn set_header(&mut self, items: &[String]) {
        (**self).set_header(items)
    }
    fn header(&self) -> &[String] {
        (**self).header()
    }
    fn write_row(&mut self, items: &[&str]) -> Result<(), VCFUtilsError> {
        (**self).write_row(items)
    }
    fn write_header(&mut self) -> Result<(), VCFUtilsError> {
        (**self).write_header()
    }
    fn write_row_bytes(&mut self, items: &[&[u8]]) -> Result<(), VCFUtilsError> {
        (**self).write_row_bytes(items)
    }
    fn write_dict(&mut self, items: &HashMap<&str, &str>) -> Result<(), VCFUtilsError> {
        (**self).write_dict(items)
    }
    fn write_dict_bytes(&mut self, items: &HashMap<&[u8], &[u8]>) -> Result<(), VCFUtilsError> {
        (**self).write_dict_bytes(items)
    }
}

#[derive(Debug)]
pub struct TSVWriter<W: Write> {
    writer: W,
    header: Vec<String>,
}

impl<W: Write> TSVWriter<W> {
    pub fn new(writer: W) -> Self {
        TSVWriter {
            writer,
            header: Vec::new(),
        }
    }
}

impl<W: Write> TableWriter for TSVWriter<W> {
    fn set_header(&mut self, items: &[String]) {
        self.header.clear();
        self.header.extend_from_slice(items);
    }

    fn header(&self) -> &[String] {
        &self.header
    }

    fn write_row(&mut self, items: &[&str]) -> Result<(), VCFUtilsError> {
        for (i, data) in items.iter().enumerate() {
            if i != 0 {
                self.writer.write_all(b"\t")?;
            }
            write!(self.writer, "{}", data)?;
        }
        writeln!(self.writer)?;
        Ok(())
    }
}

#[derive(Debug)]
pub struct CSVWriter<W: Write> {
    writer: csv::Writer<W>,
    header: Vec<String>,
}

impl<W: Write> CSVWriter<W> {
    pub fn new(writer: W) -> Self {
        CSVWriter {
            writer: csv::Writer::from_writer(writer),
            header: Vec::new(),
        }
    }
}

impl<W: Write> TableWriter for CSVWriter<W> {
    fn set_header(&mut self, items: &[String]) {
        self.header.clear();
        self.header.extend_from_slice(items);
    }

    fn header(&self) -> &[String] {
        &self.header
    }

    fn write_row(&mut self, items: &[&str]) -> Result<(), VCFUtilsError> {
        self.writer.write_record(items)?;
        Ok(())
    }
}

pub struct XlsxSheetWriter<'a, 'b> {
    writer: &'a mut xlsxwriter::Worksheet<'b>,
    header: Vec<String>,
    header_comment: Vec<String>,
    data_type: Vec<XlsxDataType>,
    current_row: u32,
    number_format: xlsxwriter::Format,
}

impl<'a, 'b> XlsxSheetWriter<'a, 'b> {
    pub fn new(workhseet: &'a mut xlsxwriter::Worksheet<'b>) -> Self {
        let mut number_format = xlsxwriter::Format::new();
        number_format.set_num_format("0.00");

        XlsxSheetWriter {
            writer: workhseet,
            header: Vec::new(),
            header_comment: Vec::new(),
            data_type: Vec::new(),
            current_row: 0,
            number_format,
        }
    }

    pub fn set_data_type(&mut self, data_type: &[XlsxDataType]) {
        self.data_type.clear();
        self.data_type.extend_from_slice(data_type);
    }

    pub fn set_header_comment(&mut self, items: &[String]) {
        self.header_comment.clear();
        self.header_comment.extend_from_slice(items);
    }

    pub fn set_column_sizes(&mut self, widths: &[f64]) -> Result<(), VCFUtilsError> {
        for (column_index, one_width) in widths.iter().enumerate() {
            self.writer
                .set_column(column_index as u16, column_index as u16, *one_width, None)?;
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum XlsxDataType {
    String,
    Boolean,
    Number,
    Formula,
}

impl<'a, 'b> TableWriter for XlsxSheetWriter<'a, 'b> {
    fn set_header(&mut self, items: &[String]) {
        self.header.clear();
        self.header.extend_from_slice(items);
    }
    fn header(&self) -> &[String] {
        &self.header
    }

    fn write_header(&mut self) -> Result<(), VCFUtilsError> {
        let header: Vec<_> = self.header().iter().map(|x| x.to_string()).collect();
        for (i, column) in header.iter().enumerate() {
            self.writer
                .write_string(self.current_row, i as u16, column, None)?;
            if let Some(comment) = self.header_comment.get(i) {
                if comment != "" {
                    self.writer
                        .write_comment(self.current_row, i as u16, &comment)?;
                }
            }
        }
        self.current_row += 1;
        Ok(())
    }

    fn write_row(&mut self, items: &[&str]) -> Result<(), VCFUtilsError> {
        for (i, column) in items.iter().enumerate() {
            if column.is_empty() {
                self.writer.write_blank(self.current_row, i as u16, None)?;
            } else {
                match self
                    .data_type
                    .get(i)
                    .copied()
                    .unwrap_or(XlsxDataType::String)
                {
                    XlsxDataType::String => {
                        if column.len() > 32766 {
                            self.writer.write_string(
                                self.current_row,
                                i as u16,
                                &format!("{}...", &column[0..32763]),
                                None,
                            )?;
                        } else {
                            self.writer
                                .write_string(self.current_row, i as u16, column, None)?;
                        }
                    }
                    XlsxDataType::Number => {
                        if let Ok(f) = column.parse() {
                            self.writer
                                .write_number(self.current_row, i as u16, f, None)?;
                        } else {
                            self.writer
                                .write_string(self.current_row, i as u16, column, None)?;
                        }
                    }
                    XlsxDataType::Boolean => {
                        self.writer.write_boolean(
                            self.current_row,
                            i as u16,
                            *column == "TRUE" || *column == "True" || *column == "true",
                            None,
                        )?;
                    }
                    XlsxDataType::Formula => {
                        self.writer.write_formula(
                            self.current_row,
                            i as u16,
                            &column,
                            Some(&self.number_format),
                        )?;
                    }
                }
            }
        }
        self.current_row += 1;
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_table_writer() -> Result<(), VCFUtilsError> {
        let mut write_buf: Vec<u8> = Vec::new();
        let mut tsv_writer = TSVWriter::new(&mut write_buf);
        tsv_writer.set_header(&["Col1".to_string(), "Col2".to_string(), "Col3".to_string()]);
        tsv_writer.write_header()?;
        tsv_writer.write_row(&["val1", "val2", "val3"])?;
        tsv_writer.write_row_bytes(&[b"bin1", b"bin2", b"bin3"])?;

        let str_map: HashMap<&str, &str> =
            vec![("Col1", "hash1"), ("Col2", "hash2"), ("Col3", "hash3")]
                .into_iter()
                .collect();
        tsv_writer.write_dict(&str_map)?;

        let bytes_map: HashMap<&[u8], &[u8]> = vec![
            (&b"Col1"[..], &b"bin_hash1"[..]),
            (&b"Col2"[..], &b"bin_hash2"[..]),
            (&b"Col3"[..], &b"bin_hash3"[..]),
        ]
        .into_iter()
        .collect();
        tsv_writer.write_dict_bytes(&bytes_map)?;

        let str_map2: HashMap<&str, &str> =
            vec![("Col1", "hash1"), ("Col3", "hash3"), ("Col4", "hash4")]
                .into_iter()
                .collect();
        tsv_writer.write_dict(&str_map2)?;

        let expected_bytes = br#"Col1	Col2	Col3
val1	val2	val3
bin1	bin2	bin3
hash1	hash2	hash3
bin_hash1	bin_hash2	bin_hash3
hash1		hash3
"#;
        assert_eq!(&expected_bytes[..], &write_buf[..]);

        Ok(())
    }

    #[test]
    fn test_csv_writer() -> Result<(), VCFUtilsError> {
        let mut write_buf: Vec<u8> = Vec::new();
        {
            let mut csv_writer = CSVWriter::new(&mut write_buf);
            csv_writer.set_header(&["Col1".to_string(), "Col2".to_string(), "Col3".to_string()]);
            csv_writer.write_header()?;
            csv_writer.write_row(&["val1", "val2,val", "val3"])?;
            csv_writer.write_row_bytes(&[b"bin1", b"bin2", b"bin3,bin4"])?;
        }

        let expected_bytes = br#"Col1,Col2,Col3
val1,"val2,val",val3
bin1,bin2,"bin3,bin4"
"#;
        assert_eq!(&expected_bytes[..], &write_buf[..]);

        Ok(())
    }
}
