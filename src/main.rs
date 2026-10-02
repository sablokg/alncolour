mod args;
use crate::args::CommandParse;
use crate::args::Commands;
use clap::Parser;
use figlet_rs::FIGfont;
mod alignmenthtml;
mod alignsubseq;

/*
Gaurav Sablok
gsablok@proton.me
*/

#[tokio::main]
async fn main() {
    let fontgenerate = FIGfont::standard().unwrap();
    let repgenerate = fontgenerate.convert("alncolor");
    println!("{}", repgenerate.unwrap());
    let args = CommandParse::parse();
    match &args.command {
        Commands::AlignmentHtml {
            inputfasta,
            outputpath,
            thread,
        } => {
            let value = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .expect("file failed");
            value.install(|| {
                let _ = alignmenthtml::alignmenthtml(inputfasta, outputpath);
                println!("The command has finished");
            });
        }
        Commands::SubAlignmentHtml {
            inputfasta,
            start,
            end,
            prefix,
            thread,
        } => {
            let value = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .expect("file failed");
            value.install(|| {
                let _ = alignsubseq::alignsubseq(inputfasta, start, end, prefix);
                println!("The command has finished");
            });
        }
    }
}
