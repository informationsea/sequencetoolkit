use super::utils::RegionMerger;
use super::BedRegion;
use std::collections::HashMap;
use std::io;

pub struct BedMerger {
    regions: HashMap<Vec<u8>, RegionMerger>,
}

impl BedMerger {
    pub fn new() -> BedMerger {
        BedMerger {
            regions: HashMap::new(),
        }
    }

    pub fn add(&mut self, region: &BedRegion) {
        if !self.regions.contains_key(&region.chromosome) {
            self.regions
                .insert(region.chromosome.clone(), RegionMerger::new());
        }
        self.regions
            .get_mut(&region.chromosome)
            .unwrap()
            .add_interval(region.start..region.end);
    }

    pub fn export_bed(&self, writer: &mut impl io::Write) -> io::Result<()> {
        let mut chromosomes: Vec<&Vec<u8>> = self.regions.keys().collect();
        chromosomes.sort();
        for one_chrom in chromosomes {
            let mut regions: Vec<_> = self
                .regions
                .get(one_chrom)
                .unwrap()
                .regions()
                .iter()
                .collect();
            regions.sort_by_cached_key(|x| x.start);
            for one_region in regions {
                writer.write_all(one_chrom)?;
                writer.write_all(b"\t")?;
                writeln!(writer, "{}\t{}", one_region.start, one_region.end)?;
            }
        }

        Ok(())
    }
}

impl Default for BedMerger {
    fn default() -> Self {
        BedMerger::new()
    }
}
