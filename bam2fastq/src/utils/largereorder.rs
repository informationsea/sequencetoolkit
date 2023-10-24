use std::fs;
use std::io::BufReader;
use std::marker::PhantomData;
use std::path::{Path, PathBuf};
use temp_file::{TempFile, TempFileBuilder};

pub struct LargeReorder<T>
where
    T: serde::Serialize + serde::de::DeserializeOwned + Clone + Ord + Eq + Sized,
{
    _phantom: PhantomData<T>,
    current_data: Vec<T>,
    maximum_data_count: usize,
    temp_dir: Option<PathBuf>,
    prefix: String,
    temp_files: Vec<TempFile>,
}

impl<T> LargeReorder<T>
where
    T: serde::Serialize + serde::de::DeserializeOwned + Clone + Ord + Eq + Sized,
{
    pub fn new(maximum_data_count: usize) -> Self {
        LargeReorder {
            _phantom: PhantomData::default(),
            current_data: Vec::new(),
            maximum_data_count,
            temp_dir: None,
            prefix: "tmp.".to_string(),
            temp_files: Vec::new(),
        }
    }

    pub fn with_temp_dir<P: AsRef<Path>>(
        maximum_data_count: usize,
        temp_dir: Option<P>,
        prefix: Option<&str>,
    ) -> Self {
        LargeReorder {
            _phantom: PhantomData::default(),
            current_data: Vec::new(),
            maximum_data_count,
            temp_dir: temp_dir.map(|x| x.as_ref().to_path_buf()),
            prefix: prefix.unwrap_or("tmp.").to_string(),
            temp_files: Vec::new(),
        }
    }

    pub fn add(&mut self, value: T) -> std::io::Result<()> {
        self.current_data.push(value);
        if self.current_data.len() >= self.maximum_data_count {
            self.write_records()?;
        }

        Ok(())
    }

    fn write_records(&mut self) -> std::io::Result<()> {
        self.current_data.sort();

        let mut temp_file_builder = TempFileBuilder::new().prefix(self.prefix.as_str());
        if let Some(temp_dir) = self.temp_dir.as_ref() {
            temp_file_builder = temp_file_builder.in_dir(&temp_dir);
        }
        let temp_file = temp_file_builder.build()?;
        let mut writer = flate2::write::GzEncoder::new(
            fs::File::create(temp_file.path())?,
            flate2::Compression::fast(),
        );

        for one in self.current_data.iter() {
            rmp_serde::encode::write(&mut writer, one)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
        }

        self.temp_files.push(temp_file);
        self.current_data.clear();
        Ok(())
    }

    pub fn into_iter(mut self) -> std::io::Result<LargeReorderIterator<T>> {
        if !self.current_data.is_empty() {
            self.write_records()?;
        }

        let mut temp_file_readers = Vec::new();
        for one in self.temp_files.iter() {
            temp_file_readers.push(BufReader::new(flate2::read::MultiGzDecoder::new(
                fs::File::open(one.path())?,
            )));
        }

        let mut current_data = Vec::new();
        for one in temp_file_readers.iter_mut() {
            let v: T = rmp_serde::from_read(one)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
            current_data.push(Some(v));
        }

        Ok(LargeReorderIterator {
            _phantom: self._phantom,
            current_data,
            temp_file_readers,
            _temp_files: self.temp_files,
        })
    }
}

pub struct LargeReorderIterator<T>
where
    T: serde::Serialize + serde::de::DeserializeOwned + Clone + Ord + Eq + Sized,
{
    _phantom: PhantomData<T>,
    _temp_files: Vec<TempFile>,
    temp_file_readers: Vec<BufReader<flate2::read::MultiGzDecoder<fs::File>>>,
    current_data: Vec<Option<T>>,
}

impl<T> LargeReorderIterator<T>
where
    T: serde::Serialize + serde::de::DeserializeOwned + Clone + Ord + Eq + Sized,
{
    fn minimum_index(&self) -> Option<usize> {
        let mut index = None;
        let mut data = None;
        for (i, one) in self.current_data.iter().enumerate() {
            if let Some(v) = one.as_ref() {
                if let Some(compare) = data {
                    if v < compare {
                        data = Some(v);
                        index = Some(i);
                    }
                } else {
                    data = Some(v);
                    index = Some(i);
                }
            }
        }
        index
    }
}

impl<T> Iterator for LargeReorderIterator<T>
where
    T: serde::Serialize + serde::de::DeserializeOwned + Clone + Ord + Eq + Sized,
{
    type Item = std::io::Result<T>;

    fn next(&mut self) -> Option<Self::Item> {
        let mi = self.minimum_index();
        if let Some(mi) = mi {
            let v = self.current_data[mi].take();

            let n: Result<T, rmp_serde::decode::Error> =
                rmp_serde::from_read(&mut self.temp_file_readers[mi]);

            match n {
                Ok(n) => {
                    self.current_data[mi] = Some(n);
                    v.map(|x| Ok(x))
                }
                Err(rmp_serde::decode::Error::InvalidMarkerRead(e))
                | Err(rmp_serde::decode::Error::InvalidDataRead(e)) => {
                    if e.kind() == std::io::ErrorKind::UnexpectedEof {
                        v.map(|x| Ok(x))
                    } else {
                        Some(Err(e))
                    }
                }
                Err(e) => Some(Err(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    e.to_string(),
                ))),
            }
        } else {
            None
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_large_reorder() -> anyhow::Result<()> {
        let mut reorder = LargeReorder::with_temp_dir(2, Some("../target"), Some("tmp."));
        for i in 0..20 {
            reorder.add(100 - i)?;
        }

        let mut result = Vec::new();
        for one in reorder.into_iter()? {
            let one = one?;
            result.push(one);
        }

        let expected: Vec<_> = (81..=100).collect();
        assert_eq!(result, expected);

        Ok(())
    }
}
