from transformers import AutoModelForCausalLM, AutoTokenizer


model_checkpoint = "as-cle-bert/saccharomyces-pythia-v1"
model = AutoModelForCausalLM.from_pretrained(model_checkpoint)
tokenizer = AutoTokenizer.from_pretrained(model_checkpoint)
