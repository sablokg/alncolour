use clap::{Parser, Subcommand};
#[derive(Debug, Parser)]
#[command(
    name = "alncolour",
    version = "1.0",
    about = " sam bam genome alignment color coded
       ************************************************
       Author Gaurav Sablok,
       Email: gsablok@proton.me
      ************************************************"
)]
pub struct CommandParse {
    /// subcommands for the specific actions
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// alignment to html
    AlignmentHtml {
        /// input path to the fasta
        inputfasta: String,
        /// output path to the html
        outputpath: String,
        /// threads for the analysis
        thread: String,
    },
    /// sub alignment to html
    SubAlignmentHtml {
        /// input path to the fasta
        inputfasta: String,
        /// start of the alignment
        start: String,
        /// end of the alignment
        end: String,
        /// prefix for the alignment
        prefix: String,
        /// threads for the analysis
        thread: String,
    },
}
