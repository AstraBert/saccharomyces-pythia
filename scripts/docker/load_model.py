from transformers import AutoModelForCausalLM, AutoTokenizer


model_checkpoint = "as-cle-bert/saccharomyces-pythia"
model = AutoModelForCausalLM.from_pretrained(model_checkpoint)
tokenizer = AutoTokenizer.from_pretrained(model_checkpoint)