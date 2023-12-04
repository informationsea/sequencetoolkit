use std::fs::File;
use std::io::Read;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

#[derive(Debug, Clone)]
pub struct Progress {
    current: Arc<AtomicUsize>,
    total: usize,
}

impl Progress {
    pub fn current(&self) -> usize {
        self.current.load(Ordering::Acquire)
    }

    pub fn total(&self) -> usize {
        self.total
    }

    pub fn progress(&self) -> f64 {
        self.current() as f64 / self.total() as f64
    }
}

pub struct ProgressReader<R: Read> {
    inner: R,
    current: Arc<AtomicUsize>,
    total: usize,
}

impl<R: Read> Read for ProgressReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let read = self.inner.read(buf)?;
        self.current.fetch_add(read, Ordering::Acquire);
        Ok(read)
    }
}

impl<R: Read> ProgressReader<R> {
    pub fn new(inner: R, total: usize) -> Self {
        ProgressReader {
            inner,
            current: Arc::new(0.into()),
            total: total,
        }
    }

    pub fn current(&self) -> usize {
        self.current.load(Ordering::Acquire)
    }

    pub fn progress(&self) -> Progress {
        Progress {
            current: self.current.clone(),
            total: self.total,
        }
    }

    pub fn total(&self) -> usize {
        self.total
    }
}

impl ProgressReader<File> {
    pub fn from_file(file: File) -> std::io::Result<Self> {
        let total = file.metadata()?.len() as usize;
        Ok(ProgressReader::new(file, total))
    }
}
